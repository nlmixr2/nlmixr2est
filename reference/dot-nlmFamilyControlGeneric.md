# Shared control setup for the nlm-family estimation methods

Shared control setup for the nlm-family estimation methods

## Usage

``` r
.nlmFamilyControlGeneric(env, controlFn, controlClass)
```

## Arguments

- env:

  dispatch environment (provides \`ui\` and \`control\`)

- controlFn:

  the method's \`\*Control()\` constructor (e.g. \`nlmControl\`)

- controlClass:

  the control object's S3 class (e.g. \`"nlmControl"\`)

## Value

Nothing; assigns the resolved control onto \`env\$ui\`

## Author

Matthew L. Fidler
