# Generic for nlmixr2 estimation methods

Generic for nlmixr2 estimation methods

## Usage

``` r
# S3 method for class 'agq'
nlmixr2Est(env, ...)

# S3 method for class 'bobyqa'
nlmixr2Est(env, ...)

# S3 method for class 'emvi'
nlmixr2Est(env, ...)

# S3 method for class 'fbvi'
nlmixr2Est(env, ...)

# S3 method for class 'fo'
nlmixr2Est(env, ...)

# S3 method for class 'foce'
nlmixr2Est(env, ...)

# S3 method for class 'focei'
nlmixr2Est(env, ...)

# S3 method for class 'output'
nlmixr2Est(env, ...)

# S3 method for class 'foceif'
nlmixr2Est(env, ...)

# S3 method for class 'focef'
nlmixr2Est(env, ...)

# S3 method for class 'focepf'
nlmixr2Est(env, ...)

# S3 method for class 'mfoceif'
nlmixr2Est(env, ...)

# S3 method for class 'mfocef'
nlmixr2Est(env, ...)

# S3 method for class 'mfocepf'
nlmixr2Est(env, ...)

# S3 method for class 'ifoceif'
nlmixr2Est(env, ...)

# S3 method for class 'ifocef'
nlmixr2Est(env, ...)

# S3 method for class 'ifocepf'
nlmixr2Est(env, ...)

# S3 method for class 'agqf'
nlmixr2Est(env, ...)

# S3 method for class 'magqf'
nlmixr2Est(env, ...)

# S3 method for class 'iagqf'
nlmixr2Est(env, ...)

# S3 method for class 'focep'
nlmixr2Est(env, ...)

# S3 method for class 'foi'
nlmixr2Est(env, ...)

# S3 method for class 'iagq'
nlmixr2Est(env, ...)

# S3 method for class 'ifoce'
nlmixr2Est(env, ...)

# S3 method for class 'ifocei'
nlmixr2Est(env, ...)

# S3 method for class 'ifocep'
nlmixr2Est(env, ...)

# S3 method for class 'ilaplace'
nlmixr2Est(env, ...)

# S3 method for class 'imp'
nlmixr2Est(env, ...)

# S3 method for class 'impmap'
nlmixr2Est(env, ...)

# S3 method for class 'laplace'
nlmixr2Est(env, ...)

# S3 method for class 'lbfgsb3c'
nlmixr2Est(env, ...)

# S3 method for class 'magq'
nlmixr2Est(env, ...)

# S3 method for class 'mfoce'
nlmixr2Est(env, ...)

# S3 method for class 'mfocei'
nlmixr2Est(env, ...)

# S3 method for class 'mfocep'
nlmixr2Est(env, ...)

# S3 method for class 'mlaplace'
nlmixr2Est(env, ...)

# S3 method for class 'n1qn1'
nlmixr2Est(env, ...)

# S3 method for class 'newuoa'
nlmixr2Est(env, ...)

# S3 method for class 'nlm'
nlmixr2Est(env, ...)

# S3 method for class 'nlme'
nlmixr2Est(env, ...)

# S3 method for class 'nlminb'
nlmixr2Est(env, ...)

nlmixr2Est(env, ...)

# Default S3 method
nlmixr2Est(env, ...)

# S3 method for class 'nls'
nlmixr2Est(env, ...)

# S3 method for class 'npag'
nlmixr2Est(env, ...)

# S3 method for class 'mnpag'
nlmixr2Est(env, ...)

# S3 method for class 'inpag'
nlmixr2Est(env, ...)

# S3 method for class 'npb'
nlmixr2Est(env, ...)

# S3 method for class 'mnpb'
nlmixr2Est(env, ...)

# S3 method for class 'inpb'
nlmixr2Est(env, ...)

# S3 method for class 'optim'
nlmixr2Est(env, ...)

# S3 method for class 'neldermead'
nlmixr2Est(env, ...)

# S3 method for class 'bfgs'
nlmixr2Est(env, ...)

# S3 method for class 'cg'
nlmixr2Est(env, ...)

# S3 method for class 'lbfgsb'
nlmixr2Est(env, ...)

# S3 method for class 'sann'
nlmixr2Est(env, ...)

# S3 method for class 'brent'
nlmixr2Est(env, ...)

# S3 method for class 'posthoc'
nlmixr2Est(env, ...)

# S3 method for class 'qrpem'
nlmixr2Est(env, ...)

# S3 method for class 'rxSolve'
nlmixr2Est(env, ...)

# S3 method for class 'simulate'
nlmixr2Est(env, ...)

# S3 method for class 'simulation'
nlmixr2Est(env, ...)

# S3 method for class 'predict'
nlmixr2Est(env, ...)

# S3 method for class 'saem'
nlmixr2Est(env, ...)

# S3 method for class 'uobyqa'
nlmixr2Est(env, ...)

# S3 method for class 'vae'
nlmixr2Est(env, ...)
```

## Arguments

- env:

  Environment for the nlmixr2 estimation routines.

  This needs to have:

  \- rxode2 ui object in \`\$ui\`

  \- data to fit in the estimation routine in \`\$data\`

  \- control for the estimation routine's control options in \`\$ui\`

- ...:

  Other arguments provided to \`nlmixr2Est()\` provided for flexibility
  but not currently used inside nlmixr

## Value

nlmixr2 fit object

## Details

This is a S3 generic that allows others to use the nlmixr2 environment
to do their own estimation routines

A prior distribution given in the \`ini()\` block is never silently
ignored. Before dispatching, \`nlmixr2Est()\` refuses any prior the
method cannot use, so a method that does not handle priors gets that for
free and a model carrying one fails with an explanation instead of being
fit to something other than what it says.

A method says what it supports with an attribute on itself:

“\` attr(nlmixr2Est.myMethod, "nlmixr2Priors") \<- "general" “\`

\- \`"none"\`, which is also what an absent attribute means – the method
cannot use priors at all - \`"theta"\` – priors on population parameters
only (\`dnorm()\`, \`dcauchy()\`, \`stdNormal()\` and the
\`multiNormal()\` family among thetas); anything referencing an omega
element is refused. No FOCEi-family method uses this level – they
declare \`"general"\`. - \`"general"\` – everything the shared kernel's
\`"general"\` method covers: the above, plus a normal prior directly on
an omega element and a textbook inverse-Wishart on an omega block. This
is what \`est="focei"\` and the rest of its family declare. -
\`"nwpri"\` – normal priors and degrees of freedom on an omega block
(\`invWishart(4)\`), evaluated the way a NONMEM NWPRI model works; a
normal prior on the omega values themselves is refused - \`"tnpri"\` –
normal priors including directly on omega elements (Monolix's/NONMEM's
own-estimation joint-normal assumption); \`dcauchy()\` and
\`invWishart()\` are refused - \`"all"\` – the method handles
everything, so nothing is checked

## Author

Matthew Fidler
