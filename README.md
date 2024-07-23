# CharacterVarieties.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://baileywhitbread.github.io/CharacterVarieties.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://baileywhitbread.github.io/CharacterVarieties.jl/dev/)
[![Build Status](https://github.com/baileywhitbread/CharacterVarieties.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/baileywhitbread/CharacterVarieties.jl/actions/workflows/CI.yml?query=branch%3Amaster)


This package computes $E$-polynomials of character varieties associated to general reductive groups in Julia using Jean Michel's package [Chevie](https://github.com/jmichel7/Chevie.jl). 

This works on Julia v1.10.2 using the Chevie version dated April 12th 2024. 

A supporting text is my [masters thesis](https://baileywhitbread.com/files/24_mphil_thesis.pdf).

## Getting started
Download and install [Julia](https://julialang.org/downloads/). In the REPL (Julia's interactive command-line), copy-paste and run the below:

```julia
using Pkg; Pkg.add(url="https://github.com/baileywhitbread/CharacterVarieties.jl")
```

This will install the `CharacterVarieties.jl` package. To load the package, copy-paste and run the below:

```julia
using CharacterVarieties
```




## Background
### Counting points
We access important cohomological information about varieties (defined by polynomials in $\mathbb{Z}[t]$) by counting points over finite fields (c.f., [HRV](https://link.springer.com/article/10.1007/s00222-008-0142-x)). If $\mathbf{A}$ is such a variety and there exists a polynomial $p_\mathbf{A}$ such that $`\#\mathbf{A}(\mathbf{F}_{q^m})=p_\mathbf{A}(q^m)`$ then we call $p_\mathbf{A}$ the $E$-polynomial of $\mathbf{A}$ and write $`E(\mathbf{A};q):=p_\mathbf{A}(q)=\#\mathbf{A}(\mathbf{F}_q)`$. Note it is not sufficient to just check $`\#\mathbf{A}(\mathbf{F}_{q})`$ is a polynomial in $q$; such polynomial must be stable under base change to finite extensions.

### Multiplicative character varieties
Let $G$ be a connected split reductive group over $\mathbb{F}_q$ with connected centre $Z$ and split maximal torus $T$. Fix integers $g\geq 0$ and $n\geq 1$, and select a strongly regular elements $S_1,\ldots,S_n$ in $T$ that are 'generic.' Let $C_1,\ldots,C_n$ be their conjugacy classes. 

The multiplicative character variety is the GIT quotient
```math
\mathbf{X} = \bigg\{(A_1,B_1,\ldots,A_g,B_g,Y_1,\ldots,Y_n)\in G^{2g}\times C_1\times \cdots\times C_n\ \bigg|\ [A_1,B_1]\cdots[A_g,B_g]Y_1\cdots Y_n = 1\bigg\}\bigg/\!\!\!\!\!\bigg/G
```
where the action is simultaneous conjugation. 

This package computes the $E$-polynomial $E(\mathbf{X};q)$ via our formula for $`\#\mathbf{X}(\mathbb{F}_q)`$.  

### Additive character varieties
Let $\mathfrak{g}$ be the Lie algebra of $G$ and let $\mathfrak{t}$ be the Lie algebra of $T$. Select regular semisimple elements $s_1,\ldots,s_n$ in $\mathfrak{t}$ that are 'generic.' Let $O_1,\ldots,O_n$ be their adjoint orbits. 

The additive character variety is the GIT quotient
```math
\mathbf{Y} = \bigg\{(A_1,B_1,\ldots,A_g,B_g,Y_1,\ldots,Y_n)\in \mathfrak{g}^{2g}\times O_1\times \cdots\times O_n\ \bigg|\ [A_1,B_1]+\cdots+[A_g,B_g] + Y_1+ \cdots + Y_n = 0\bigg\}\bigg/\!\!\!\!\!\bigg/G
```
where the action is simultaneous conjugation (i.e., the adjoint action). 

This package computes the $E$-polynomial $E(\mathbf{Y};q)$ via our formula for $`\#\mathbf{Y}(\mathbb{F}_q)`$.  


## Calculating E-polynomials
We will use the group $G=G_2$ (i.e., the semisimple group of adjoint type $G_2$) as an example. To select this group, we use the command `G=rootdatum(:G2)`. One could have instead chosen, for instance, `rootdatum(:gl,2)`, `rootdatum(:so,5)`, `rootdatum(:pgl,3)`, or `rootdatum(:F4)`.

Suppose $g\geq 0$ is the genus number and $n\geq 1$ is the number of punctures. Then
- The command `EX(G,g,n)` returns the $E$-polynomial $E(\mathbf{X};q)$, and
- The command `EY(G,g,n)` returns the $E$-polynomial $E(\mathbf{Y};q)$.

For instance, 
```julia
julia> EX(G,0,3)
Frac{Pol{Rational{Int64}}}: q⁸+6q⁷+20q⁶+58q⁵+180q⁴+58q³+20q²+6q+1

julia> EY(G,0,3)
Pol{BigInt}: q⁸+6q⁷+19q⁶+45q⁵+99q⁴
```



## To do
- Clean up polynomial types their coefficient types. I used `BigInt` to verify non-negativity of coefficients of $E(\mathbf{Y};q)$ for large $G$,$g$ and $n$. This will probably require calculating polynomials more efficiently.
- Speed up $E$-polynomial calculations. The major killer is calculating $G$-types and their associated data. For instance, just calculating the $F_4$-types takes almost 10 seconds:
```julia
@time group_types(rootdatum(:F4))
  8.824810 seconds (11.60 M allocations: 787.033 MiB, 2.69% gc time, 99.45% compilation time)
```
- Add real-time calculations of $E$-polynomials on my website [baileywhitbread.com](https://www.baileywhitbread.com).




## Further directions

### Mixed Hodge polynomials
Associated to $\mathbf{X}$ is the (compactly supported) mixed Hodge polynomial
```math
H(\mathbf{X};x,y,t) = \sum_{i,j,k} h^{i,j,k} x^i y^j t^k.
```
Proving that $`\#\mathbf{X}(\mathbb{F}_q)`$ is a polynomial in $q$ implies $H(\mathbf{X};x,y,t)$ depends only on the product $xy=:q$ and $t$, and that
```math
H(\mathbf{X};q,-1) = E(\mathbf{X};q).
```
That is, computing $`\#\mathbf{X}(\mathbb{F}_q)`$ gives us a specialisation of $H(\mathbf{X};q,t)$. 

### Pure parts
Contained in the cohomology ring $H^\ast(\mathbf{X})$ is an important subring called the pure subring $H^\ast_\mathrm{pure}(\mathbf{X})$. Denote by $PH(\mathbf{X};u)$ the Poincar\'{e} polynomial of the pure subring. Then there is another important specialisation $H(\mathbf{X})\mapsto PH(\mathbf{X})$ defined by setting all terms to zero except the monomials in $qt^2$. For instance, 
```math
H(\mathbf{X};q,t)=qt^2 + q^2t^4 + q + qt + 1 \rightsquigarrow PH(\mathbf{X};q)=q+q^2+1.
```
We say a variety's cohomology is 'pure' if all terms are monomials in $qt^2$. 

### Non-negative coefficients
When $G=\mathrm{GL}_n$, it was proven that the coefficients of $`\#\mathbf{Y}(\mathbb{F}_q)`$ are non-negative by relating the additive character variety to a quiver variety (which are known to have pure cohomology due to their symplectic geometry, c.f., [HLRV](https://projecteuclid.org/journals/duke-mathematical-journal/volume-160/issue-2/Arithmetic-harmonic-analysis-on-character-and-quiver-varieties/10.1215/00127094-1444258.full)). Moreover, it was conjectured the polynomials $\mathbf{X}$ and $\mathbf{X}$ are closely related, in the sense that

```math
PH(\mathbf{X};q) = E(\mathbf{Y};q).
```

The conjecture has been proven in one narrow case: $G=\mathrm{GL}_2$, $n=1$ and $C_1$ is the conjugacy class of $`\left(\begin{smallmatrix}-1 & \\ & -1 \end{smallmatrix}\right)`$. An unproven conjectural formula is given for $G=\mathrm{GL}_n$ in [HLRV](https://projecteuclid.org/journals/duke-mathematical-journal/volume-160/issue-2/Arithmetic-harmonic-analysis-on-character-and-quiver-varieties/10.1215/00127094-1444258.full). 


### An idea
Use the specialisations $E(\mathbf{X};q)$ and $E(\mathbf{Y};q)$ to try and understand the mixed Hodge polynomial $H(\mathbf{X};q,t)$.

### Another idea
Use the `CharacterVarieties.jl` to search for $E$-polynomials of additive character varieties with negative coefficients.

### Another another idea
If $G=\mathrm{GL}_n$ then the formulas for $E(\mathbf{X};q)$ and $E(\mathbf{Y};q)$ in [HLRV](https://projecteuclid.org/journals/duke-mathematical-journal/volume-160/issue-2/Arithmetic-harmonic-analysis-on-character-and-quiver-varieties/10.1215/00127094-1444258.full) are very different to ours. It would be interesting to implement their formulas, especially given their conjectural formula for the mixed Hodge polynomial. These formulas involve symmetric functions and inner products of complete and monomial symmetric functions.  