suppressMessages(library(symengine))
S <- symengine::S; d <- function(e,v) symengine::D(e,S(v))
PHI <- function(a) sprintf("0.5*(1+erf((%s)/sqrt(2)))", a)
rho <- list(
  M3  = S(sprintf("-log(%s)", PHI("cens*(dv-f)/sqrt(r)"))),
  M4  = S(sprintf("-(log(%s - %s) - log(1 - %s))",
             PHI("cens*(dv-f)/sqrt(r)"), PHI("cens*(lim-f)/sqrt(r)"), PHI("cens*(lim-f)/sqrt(r)"))),
  M2g = S(sprintf("0.5*((dv-f)^2/r + log(r)) + log(%s)", PHI("-(lim-f)/sqrt(r)"))),
  M2l = S(sprintf("0.5*((dv-f)^2/r + log(r)) + log(%s)", PHI("(lim-f)/sqrt(r)"))))
specs2 <- list(f="f", r="r", ff=c("f","f"), fr=c("f","r"), rr=c("r","r"))
specs3 <- list(fff=c("f","f","f"), ffr=c("f","f","r"), frr=c("f","r","r"), rrr=c("r","r","r"))
idx <- c(f=0,r=1,ff=2,fr=3,rr=4,fff=5,ffr=6,frr=7,rrr=8)
partial <- function(e, sp){ for(v in sp) e<-d(e,v); e }
r2c <- function(node) {
  if (is.numeric(node)) return(formatC(node, format="f", digits=15))
  if (is.name(node)) { nm<-as.character(node); return(if(nm=="pi") "M_PI" else nm) }
  if (is.call(node)) { op <- as.character(node[[1]])
    if (length(all.vars(node))==0) return(formatC(eval(node), format="f", digits=15))
    if (op=="(") return(sprintf("(%s)", r2c(node[[2]])))
    if (op=="/" && length(node)==3) return(sprintf("(%s / _safe_zero(%s))", r2c(node[[2]]), r2c(node[[3]])))
    if (op %in% c("+","-","*") && length(node)==3) return(sprintf("(%s %s %s)", r2c(node[[2]]), op, r2c(node[[3]])))
    if (op=="-" && length(node)==2) return(sprintf("(-%s)", r2c(node[[2]])))
    if (op=="^") { b<-node[[3]]; if (is.numeric(b)&&b==2) return(sprintf("(%s*%s)", r2c(node[[2]]),r2c(node[[2]])))
      return(sprintf("R_pow(%s, %s)", r2c(node[[2]]), r2c(b))) }
    if (op %in% c("sqrt","exp","log","erf")) return(sprintf("%s(%s)", op, r2c(node[[2]])))
    stop("unhandled fn: ", op) }
  stop("unhandled") }
genBlock <- function(m, specs, indent="    ") {
  blk <- paste(vapply(names(specs), function(sp)
           sprintf("d_%s=%s", sp, as.character(partial(rho[[m]], specs[[sp]]))), character(1)), collapse="\n")
  opt <- suppressMessages({ tf<-tempfile(); sink(tf); on.exit(sink()); r<-rxode2::rxOptExpr(blk,"c"); sink(); r })
  lines <- strsplit(opt, "\n")[[1]]; lines <- lines[nzchar(trimws(lines)) & !grepl("====|→|optimiz|duplicate", lines)]
  out <- character(0)
  for (ln in lines) {
    lr <- regmatches(ln, regexec("^\\s*([A-Za-z0-9_]+)\\s*[~=]\\s*(.*)$", ln))[[1]]
    if (length(lr)<3 || !nzchar(lr[2])) next
    lhs <- lr[2]; rhs <- r2c(parse(text=lr[3])[[1]])
    if (startsWith(lhs, "d_")) out <- c(out, sprintf("%sout[%d] = %s;", indent, idx[[sub("^d_","",lhs)]], rhs))
    else out <- c(out, sprintf("%sdouble %s = %s;", indent, lhs, rhs)) }
  paste(out, collapse="\n") }
# compute ALL blocks first (rxOptExpr progress happens here, not into the file)
B <- list(
  M2l2=genBlock("M2l",specs2,"      "), M2l3=genBlock("M2l",specs3,"        "),
  M2g2=genBlock("M2g",specs2,"      "), M2g3=genBlock("M2g",specs3,"        "),
  M42 =genBlock("M4", specs2,"    "),   M43 =genBlock("M4", specs3,"      "),
  M32 =genBlock("M3", specs2,"    "),   M33 =genBlock("M3", specs3,"      "))
con <- file("censPartials_gen.c","w")
wl <- function(...) writeLines(paste0(...), con)
wl("// GENERATED -- censored rho(f,r) partials; out[0..8]=rho_{f,r,ff,fr,rr,fff,ffr,frr,rrr}; rho=-logLik")
wl("static inline void censNormalPartials(double cens, double limDv, double lim,")
wl("                                      double f, double r, int order, double* out) {")
wl("  double dv = limDv;")
wl("  int hasFin = (R_FINITE(lim) && !ISNA(lim));")
wl("  if (cens == 0.0 && hasFin) {            // M2")
wl("    if (lim < f) {"); wl(B$M2l2); wl("      if (order >= 3) {"); wl(B$M2l3); wl("      }")
wl("    } else {");        wl(B$M2g2); wl("      if (order >= 3) {"); wl(B$M2g3); wl("      }"); wl("    }")
wl("  } else if (hasFin) {                    // M4"); wl(B$M42); wl("    if (order >= 3) {"); wl(B$M43); wl("    }")
wl("  } else {                                // M3"); wl(B$M32); wl("    if (order >= 3) {"); wl(B$M33); wl("    }")
wl("  }")
wl("}")
close(con)
cat("wrote", length(readLines("censPartials_gen.c")), "lines\n")
