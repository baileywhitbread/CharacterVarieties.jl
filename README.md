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




## Background
### Counting points
We say an algebraic variety $\mathbf{A}$ defined by polynomials in $\mathbb{Z}[t]$ is polynomial count if there exists a polynomial $p_\mathbf{A}$ such that $`\#\mathbf{A}(\mathbf{F}_{q^m})=p_\mathbf{A}(q^m)`$ for all $m\geq 1$. We call $p_\mathbf{A}$ the $E$-polynomial of $\mathbf{A}$ and write $`E(\mathbf{A};q):=p_\mathbf{A}(q)=\#\mathbf{A}(\mathbf{F}_q)`$ (c.f., [HRV](https://link.springer.com/article/10.1007/s00222-008-0142-x)).

### Multiplicative and additive character varieties
Let $G$ be a connected split reductive group over $\mathbb{F}_q$ with connected centre $Z$ and split maximal torus $T$. Fix integers $g\geq 0$ and $n\geq 1$, and select a strongly regular elements $S_1,\ldots,S_n$ in $T$ that are 'generic.' Let $C_1,\ldots,C_n$ be their conjugacy classes. 

The multiplicative character variety is the GIT quotient
```math
\mathbf{X} = \bigg\{(A_1,B_1,\ldots,A_g,B_g,Y_1,\ldots,Y_n)\in G^{2g}\times C_1\times \cdots\times C_n\ \bigg|\ [A_1,B_1]\cdots[A_g,B_g]Y_1\cdots Y_n = 1\bigg\}\bigg/\!\!\!\!\!\bigg/G
```
where the action is simultaneous conjugation. 

Let $\mathfrak{g}$ be the Lie algebra of $G$ and let $\mathfrak{t}$ be the Lie algebra of $T$. Select regular semisimple elements $s_1,\ldots,s_n$ in $\mathfrak{t}$ that are 'generic.' Let $O_1,\ldots,O_n$ be their adjoint orbits. 

The additive character variety is the GIT quotient
```math
\mathbf{Y} = \bigg\{(A_1,B_1,\ldots,A_g,B_g,Y_1,\ldots,Y_n)\in \mathfrak{g}^{2g}\times O_1\times \cdots\times O_n\ \bigg|\ [A_1,B_1]+\cdots+[A_g,B_g] + Y_1+ \cdots + Y_n = 0\bigg\}\bigg/\!\!\!\!\!\bigg/G
```
where the action is simultaneous conjugation (i.e., the adjoint action). 

## Calculating E-polynomials
This package computes the $E$-polynomials $E(\mathbf{X};q)$ and $E(\mathbf{Y};q)$ via our formulas for $`\#\mathbf{X}(\mathbb{F}_q)`$ and $`\#\mathbf{Y}(\mathbb{F}_q)`$. We will use the group $G=G_2$ (i.e., the semisimple group of adjoint type $G_2$) as an example. To select this group, we use the command `G=rootdatum(:G2)`. One can instead choose `rootdatum(:gl,2)`, `rootdatum(:so,5)`, `rootdatum(:pgl,3)`, `rootdatum(:F4)`, etc.

- The command `EX(G,g,n)` returns the $E$-polynomial $E(\mathbf{X};q)$, and
- The command `EY(G,g,n)` returns the $E$-polynomial $E(\mathbf{Y};q)$.

When $g=0$, the `Int64` data type is insufficient for polynomial division. In this case, we can use the `bigint_EX` and `bigint_EY` functions:

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

### Pure parts
Contained in the cohomology ring $H^\ast(\mathbf{X})$ is an important subring called the pure subring $H^\ast_\mathrm{pure}(\mathbf{X})$. Denote by $PH(\mathbf{X};u)$ the Poincar\'{e} polynomial of the pure subring. There' an important specialisation $H(\mathbf{X})\mapsto PH(\mathbf{X})$ given by setting all terms to zero except monomials in $qt^2$:
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