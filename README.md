# CharacterVarieties.jl

This package computes $E$-polynomials of character varieties associated to general reductive groups in Julia using Jean Michel's package [Chevie](https://github.com/jmichel7/Chevie.jl). A supporting text is my [masters thesis](https://baileywhitbread.com/files/24_mphil_thesis.pdf).





## Getting started
Download and install [Julia](https://julialang.org/downloads/). In the REPL (Julia's interactive command-line), copy-paste and run the below:

```julia
using Pkg; Pkg.add(url="https://github.com/baileywhitbread/CharacterVarieties.jl")
```

This will install the `CharacterVarieties.jl` package. To load the package, copy-paste and run the below:

```julia
using CharacterVarieties
```





## Multiplicative and additive character varieties
Fix integers $g\geq 0$ and $n\geq 1$ and let $G$ be a connected split reductive group over $\mathbb{F}_q$ with connected centre $Z$ and split maximal torus $T$. 

Select strongly regular elements $S_1,\ldots,S_n$ in $T$ that are 'generic' and let $C_1,\ldots,C_n$ be their conjugacy classes.

The multiplicative character variety is the GIT quotient
```math
\mathbf{X} = \bigg\{(A_1,B_1,\ldots,A_g,B_g,Y_1,\ldots,Y_n)\in G^{2g}\times C_1\times \cdots\times C_n\ \bigg|\ [A_1,B_1]\cdots[A_g,B_g]Y_1\cdots Y_n = 1\bigg\}\bigg/\!\!\!\!\!\bigg/G
```
where the action is simultaneous conjugation. 

Let $\mathfrak{g}$ be the Lie algebra of $G$ and let $\mathfrak{t}$ be the Lie algebra of $T$. 

Select regular semisimple elements $s_1,\ldots,s_n$ in $\mathfrak{t}$ that are 'generic' and let $O_1,\ldots,O_n$ be their adjoint orbits. 

The additive character variety is the GIT quotient
```math
\mathbf{Y} = \bigg\{(A_1,B_1,\ldots,A_g,B_g,Y_1,\ldots,Y_n)\in \mathfrak{g}^{2g}\times O_1\times \cdots\times O_n\ \bigg|\ [A_1,B_1]+\cdots+[A_g,B_g] + Y_1+ \cdots + Y_n = 0\bigg\}\bigg/\!\!\!\!\!\bigg/G
```
where the action is simultaneous conjugation (i.e., the adjoint action). 

## Calculating E-polynomials
This package computes the $E$-polynomials $E(\mathbf{X};q)$ and $E(\mathbf{Y};q)$ (see [HRV](https://link.springer.com/article/10.1007/s00222-008-0142-x) for the definition of $E$-polynomials). 

This is done using our formulas for $E(\mathbf{X};q)$ and $E(\mathbf{Y};q)$.

We will use the group $G=G_2$ (i.e., the semisimple group of adjoint type $G_2$) as an example. 

The command `G=rootdatum(:G2)` selects the group. 

(One can instead choose `rootdatum(:gl,2)`, `rootdatum(:so,5)`, `rootdatum(:pgl,3)`, `rootdatum(:F4)`, etc.)

- `EX(G,g,n)` returns $E(\mathbf{X};q)$, and
- `EY(G,g,n)` returns $E(\mathbf{Y};q)$.

When $g=0$, the `Int64` data type produces rounding errors; we use `bigint_EX` and `bigint_EY` instead:

```julia
julia> EX(G,0,3)
Pol{Int64}: q⁸+6q⁷+20q⁶+58q⁵+180q⁴+58q³+20q²+6q+1

julia> bigint_EX(G,0,3)
Pol{BigInt}: q⁸+6q⁷+20q⁶+58q⁵+180q⁴+58q³+20q²+6q+1

julia> EY(G,0,3)
ERROR: cannot convert Frac(Pol(BigFloat[2.77555756156289135105907917022705078125e-17, 0.0, 3.3306690738754696212708950042724609375e-16, 5.5511151231257827021181583404541015625e-16, -99.0, -45.000000000000000610622663543836097232997417449951171875, 80.0, 39.0, 18.0, 6.0, 1.0]),Pol(BigFloat[-1.0, 0.0, 1.0])) to Pol

julia> bigint_EY(G,0,3)
Pol{BigInt}: q⁸+6q⁷+19q⁶+45q⁵+99q⁴
```






## To do
- Speed up $E$-polynomial calculations:
```julia
julia> @time group_types(rootdatum(:F4))
 38.772694 seconds (91.92 M allocations: 6.047 GiB, 5.86% gc time, 99.80% compilation time)
julia> @time EX(rootdatum(:F4),1,1)
61.184813 seconds (227.91 M allocations: 18.563 GiB, 8.07% gc time, 88.52% compilation time)
```


- There is a large, unexplained leap in computation time between $\mathrm{SO}_9$ and $\mathrm{SO}_{11}$:
```julia
julia> @time EX(rootdatum(:so,7),0,3)
  0.622608 seconds (1.15 M allocations: 88.349 MiB, 96.49% compilation time)

julia> @time EX(rootdatum(:so,9),0,3)
  0.580247 seconds (7.36 M allocations: 715.632 MiB, 12.29% gc time, 28.14% compilation time)

julia> @time EX(rootdatum(:so,11),0,3)
127.007951 seconds (2.15 G allocations: 218.391 GiB, 17.85% gc time, 0.44% compilation time)
```


- Add real-time calculations of $E$-polynomials to [baileywhitbread.com](https://www.baileywhitbread.com).




## Further directions

### Mixed Hodge polynomials
Associated to $\mathbf{X}$ is the (compactly supported) mixed Hodge polynomial
```math
H(\mathbf{X};x,y,t) = \sum_{i,j,k} h^{i,j,k} x^i y^j t^k.
```
Proving $\mathbf{X}$ is polynomial count implies $H(\mathbf{X};x,y,t)$ depends only on the product $xy=:q$ and $t$, and that
```math
H(\mathbf{X};q,-1) = E(\mathbf{X};q).
```
There's another specialisation of $H(\mathbf{X};q,t)$ given by setting all terms to zero except monomials in $u:=qt^2$:
```math
H(\mathbf{X};q,t)=qt^2 + q^2t^4 + q + qt + 1 \rightsquigarrow PH(\mathbf{X};u):=u+u^2+1.
```
It is conjectured the polynomials $H(\mathbf{X};q,t)$ and $E(\mathbf{Y};u)$ are closely related, in the sense that
```math
PH(\mathbf{X};q) = E(\mathbf{Y};q).
```
The conjecture has been proven in one narrow case: $G=\mathrm{GL}_2$, $n=1$ and $C_1$ is the conjugacy class of $`\left(\begin{smallmatrix}-1 & \\ & -1 \end{smallmatrix}\right)`$. This was possible because $H(\mathbf{X};q,t)$ is explicitly known. When $G=\mathrm{GL}_d$ and the $C_i$ are semisimple, an unproven formula for $H(\mathbf{X};q,t)$ was given in [HLRV](https://projecteuclid.org/journals/duke-mathematical-journal/volume-160/issue-2/Arithmetic-harmonic-analysis-on-character-and-quiver-varieties/10.1215/00127094-1444258.full). 


#### Idea: Use the specialisations $E(\mathbf{X};q)$ and $E(\mathbf{Y};q)$ to try and understand the mixed Hodge polynomial $H(\mathbf{X};q,t)$ for general reductive $G$.

The formulas for $E(\mathbf{X};q)$ and $E(\mathbf{Y};q)$ in [HLRV](https://projecteuclid.org/journals/duke-mathematical-journal/volume-160/issue-2/Arithmetic-harmonic-analysis-on-character-and-quiver-varieties/10.1215/00127094-1444258.full) are very different to ours. 

It would be interesting to implement their formulas, especially given their conjectural formula for the mixed Hodge polynomial, which involve symmetric functions and inner products of complete and monomial symmetric functions.  

### Non-negative coefficients
When $G=\mathrm{GL}_d$, the coefficients of $`\#\mathbf{Y}(\mathbb{F}_q)`$ are non-negative by the work of [HLRV](https://projecteuclid.org/journals/duke-mathematical-journal/volume-160/issue-2/Arithmetic-harmonic-analysis-on-character-and-quiver-varieties/10.1215/00127094-1444258.full).

#### Idea: Use the `CharacterVarieties.jl` to search for $E(\mathbf{Y};q)$ with negative coefficients.
