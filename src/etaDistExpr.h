#ifndef __ETADISTEXPR_H__
#define __ETADISTEXPR_H__
// Evaluate a declared distribution's ARGUMENT expressions in C++.
//
// dist(eta.cl) ~ dgamma(shape = 1/exp(lclrv), rate = 1/(exp(lclrv)*exp(lclm)))
// gives arguments that are arbitrary expressions over ini() thetas.  The M-step
// fits the family's NATIVE parameters and then has to map them back onto those
// thetas, which is a small inverse problem in the thetas -- and evaluating the
// expressions is the only part that needed R.
//
// That callback was the last piece of the M-step living outside C++: one
// eval() per objective evaluation, inside an R Nelder-Mead, called from the
// C++ M-step loop.  This removes it for the expressions these declarations
// actually use; anything outside the grammar below returns false from
// etaDistExprParse and the caller keeps the R route, which is always correct.
//
// Parsed ONCE into RPN, then evaluated many times by the optimizer -- so the
// per-evaluation cost is a walk over a small vector, not a re-parse.
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>

typedef enum {
  etaDistOpNum = 0, etaDistOpVar, etaDistOpAdd, etaDistOpSub, etaDistOpMul,
  etaDistOpDiv, etaDistOpPow, etaDistOpNeg, etaDistOpExp, etaDistOpLog,
  etaDistOpSqrt, etaDistOpAbs
} etaDistOp;

typedef struct {
  etaDistOp op;
  double num;   // etaDistOpNum
  int var;      // etaDistOpVar: index into the theta vector
} etaDistTok;

// ---- tokenizer + shunting-yard ---------------------------------------------
struct etaDistLex {
  const std::string &s;
  size_t i;
  etaDistLex(const std::string &str) : s(str), i(0) {}
  void ws() { while (i < s.size() && isspace((unsigned char)s[i])) ++i; }
  bool eof() { ws(); return i >= s.size(); }
  char peek() { ws(); return i < s.size() ? s[i] : '\0'; }
};

static inline int etaDistPrec(char c) {
  switch (c) {
  case '+': case '-': return 1;
  case '*': case '/': return 2;
  case '^': return 3;
  default: return 0;
  }
}

// Returns false when the expression uses anything outside the grammar, which is
// the signal to fall back rather than to guess.
static inline bool etaDistExprParse(const std::string &expr,
                                    const std::vector<std::string> &vars,
                                    std::vector<etaDistTok> &out) {
  out.clear();
  std::vector<char> ops;         // operators and '(' ; 'e','l','s','a','n' = fns/neg
  etaDistLex lx(expr);
  bool wantVal = true;           // distinguishes unary minus from binary
  while (!lx.eof()) {
    char c = lx.peek();
    if (wantVal && (c == '-' || c == '+')) {
      lx.i++;
      if (c == '-') ops.push_back('n');
      continue;
    }
    if (isdigit((unsigned char)c) || c == '.') {
      size_t j = lx.i; char *end = NULL;
      double v = strtod(lx.s.c_str() + j, &end);
      if (end == lx.s.c_str() + j) return false;
      lx.i = (size_t)(end - lx.s.c_str());
      etaDistTok t; t.op = etaDistOpNum; t.num = v; t.var = -1;
      out.push_back(t); wantVal = false; continue;
    }
    if (isalpha((unsigned char)c) || c == '.' || c == '_') {
      size_t j = lx.i;
      while (lx.i < lx.s.size() &&
             (isalnum((unsigned char)lx.s[lx.i]) || lx.s[lx.i] == '.' ||
              lx.s[lx.i] == '_')) lx.i++;
      std::string nm = lx.s.substr(j, lx.i - j);
      lx.ws();
      if (lx.i < lx.s.size() && lx.s[lx.i] == '(') {
        lx.i++;
        if (nm == "exp") ops.push_back('e');
        else if (nm == "log") ops.push_back('l');
        else if (nm == "sqrt") ops.push_back('s');
        else if (nm == "abs") ops.push_back('a');
        else return false;       // unsupported function -> fall back
        ops.push_back('(');
        wantVal = true; continue;
      }
      int vi = -1;
      for (size_t k = 0; k < vars.size(); ++k) if (vars[k] == nm) { vi = (int)k; break; }
      if (vi < 0) return false;  // unknown symbol -> fall back
      etaDistTok t; t.op = etaDistOpVar; t.num = 0.0; t.var = vi;
      out.push_back(t); wantVal = false; continue;
    }
    if (c == '(') { lx.i++; ops.push_back('('); wantVal = true; continue; }
    if (c == ')') {
      lx.i++;
      bool found = false;
      while (!ops.empty()) {
        char o = ops.back(); ops.pop_back();
        if (o == '(') { found = true; break; }
        etaDistTok t; t.num = 0.0; t.var = -1;
        switch (o) {
        case '+': t.op = etaDistOpAdd; break;  case '-': t.op = etaDistOpSub; break;
        case '*': t.op = etaDistOpMul; break;  case '/': t.op = etaDistOpDiv; break;
        case '^': t.op = etaDistOpPow; break;  case 'n': t.op = etaDistOpNeg; break;
        default: return false;
        }
        out.push_back(t);
      }
      if (!found) return false;
      // a function application closes with its own operator
      if (!ops.empty()) {
        char o = ops.back();
        if (o == 'e' || o == 'l' || o == 's' || o == 'a') {
          ops.pop_back();
          etaDistTok t; t.num = 0.0; t.var = -1;
          t.op = (o == 'e') ? etaDistOpExp : (o == 'l') ? etaDistOpLog :
                 (o == 's') ? etaDistOpSqrt : etaDistOpAbs;
          out.push_back(t);
        }
      }
      wantVal = false; continue;
    }
    if (c == '+' || c == '-' || c == '*' || c == '/' || c == '^') {
      lx.i++;
      while (!ops.empty() && ops.back() != '(' &&
             etaDistPrec(ops.back()) >= etaDistPrec(c) && c != '^') {
        char o = ops.back(); ops.pop_back();
        etaDistTok t; t.num = 0.0; t.var = -1;
        switch (o) {
        case '+': t.op = etaDistOpAdd; break;  case '-': t.op = etaDistOpSub; break;
        case '*': t.op = etaDistOpMul; break;  case '/': t.op = etaDistOpDiv; break;
        case '^': t.op = etaDistOpPow; break;  case 'n': t.op = etaDistOpNeg; break;
        default: return false;
        }
        out.push_back(t);
      }
      ops.push_back(c); wantVal = true; continue;
    }
    return false;                // anything else -> fall back
  }
  while (!ops.empty()) {
    char o = ops.back(); ops.pop_back();
    if (o == '(') return false;
    etaDistTok t; t.num = 0.0; t.var = -1;
    switch (o) {
    case '+': t.op = etaDistOpAdd; break;  case '-': t.op = etaDistOpSub; break;
    case '*': t.op = etaDistOpMul; break;  case '/': t.op = etaDistOpDiv; break;
    case '^': t.op = etaDistOpPow; break;  case 'n': t.op = etaDistOpNeg; break;
    case 'e': t.op = etaDistOpExp; break;  case 'l': t.op = etaDistOpLog; break;
    case 's': t.op = etaDistOpSqrt; break; case 'a': t.op = etaDistOpAbs; break;
    default: return false;
    }
    out.push_back(t);
  }
  return !out.empty();
}

// NaN on a malformed stack rather than reading past it.
static inline double etaDistExprEval(const std::vector<etaDistTok> &rpn,
                                     const double *vals, int nvals) {
  double st[64]; int sp = 0;
  for (size_t i = 0; i < rpn.size(); ++i) {
    const etaDistTok &t = rpn[i];
    if (t.op == etaDistOpNum) { if (sp >= 64) return NAN; st[sp++] = t.num; continue; }
    if (t.op == etaDistOpVar) {
      if (sp >= 64 || t.var < 0 || t.var >= nvals) return NAN;
      st[sp++] = vals[t.var]; continue;
    }
    if (t.op == etaDistOpNeg || t.op == etaDistOpExp || t.op == etaDistOpLog ||
        t.op == etaDistOpSqrt || t.op == etaDistOpAbs) {
      if (sp < 1) return NAN;
      double a = st[sp-1];
      switch (t.op) {
      case etaDistOpNeg:  st[sp-1] = -a; break;
      case etaDistOpExp:  st[sp-1] = std::exp(a); break;
      case etaDistOpLog:  st[sp-1] = std::log(a); break;
      case etaDistOpSqrt: st[sp-1] = std::sqrt(a); break;
      default:            st[sp-1] = std::fabs(a); break;
      }
      continue;
    }
    if (sp < 2) return NAN;
    double b = st[--sp], a = st[sp-1];
    switch (t.op) {
    case etaDistOpAdd: st[sp-1] = a + b; break;
    case etaDistOpSub: st[sp-1] = a - b; break;
    case etaDistOpMul: st[sp-1] = a * b; break;
    case etaDistOpDiv: st[sp-1] = a / b; break;
    case etaDistOpPow: st[sp-1] = std::pow(a, b); break;
    default: return NAN;
    }
  }
  return (sp == 1) ? st[0] : NAN;
}

#endif // __ETADISTEXPR_H__
