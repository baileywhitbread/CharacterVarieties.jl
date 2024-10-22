# CharacterVarieties.jl

This package computes $E$-polynomials of character varieties associated to general reductive groups in Julia. This was written for the paper [arXiv:2409.04735](https://arxiv.org/abs/2409.04735). We heavily rely on Jean Michel's computer algebra system [Chevie](https://github.com/jmichel7/Chevie.jl). 





## Installing
Download and install [Julia](https://julialang.org/downloads/). In the REPL (Julia's interactive command-line), copy-paste and run the below:

```julia
using Pkg; Pkg.add(url="https://github.com/baileywhitbread/CharacterVarieties.jl")
```

This installs `CharacterVarieties.jl` and its dependencies. To load the package, copy-paste and run the below:

```julia
using CharacterVarieties
```





## Multiplicative and additive character varieties
Fix integers $g\geq 0$ and $n\geq 1$ and let $G$ be a connected split reductive group over $\mathbb{F}_q$ with connected centre $Z$ and split maximal torus $T$. Select strongly regular elements $S_1,\ldots,S_n$ in $T$ that are 'generic' and let $C_1,\ldots,C_n$ be their conjugacy classes. The multiplicative character variety is the GIT quotient
```math
\mathbf{X} = \bigg\{(A_1,B_1,\ldots,A_g,B_g,Y_1,\ldots,Y_n)\in G^{2g}\times C_1\times \cdots\times C_n\ \bigg|\ [A_1,B_1]\cdots[A_g,B_g]Y_1\cdots Y_n = 1\bigg\}\bigg/\!\!\!\!\!\bigg/G
```
where the action is simultaneous conjugation. This is an affine scheme of finite type over the finite field of size $q$. Let $\mathfrak{g}$ be the Lie algebra of $G$ and let $\mathfrak{t}$ be the Lie algebra of $T$. Select regular semisimple elements $s_1,\ldots,s_n$ in $\mathfrak{t}$ that are 'generic' and let $O_1,\ldots,O_n$ be their adjoint orbits. The additive character variety is the GIT quotient
```math
\mathbf{Y} = \bigg\{(A_1,B_1,\ldots,A_g,B_g,Y_1,\ldots,Y_n)\in \mathfrak{g}^{2g}\times O_1\times \cdots\times O_n\ \bigg|\ [A_1,B_1]+\cdots+[A_g,B_g] + Y_1+ \cdots + Y_n = 0\bigg\}\bigg/\!\!\!\!\!\bigg/G
```
where the action is simultaneous conjugation (i.e., the adjoint action). This is an affine scheme of finite type over the finite field of size $q$.

## Calculating E-polynomials
This package computes the $E$-polynomials $E(\mathbf{X};q)$ and $E(\mathbf{Y};q)$ (see [LRV](https://aif.centre-mersenne.org/articles/10.5802/aif.3540/) for the definition of $E$-polynomials of varieties over finite fields or [HRV](https://link.springer.com/article/10.1007/s00222-008-0142-x) for the complex analogue). This is done using our formulas for $E(\mathbf{X};q)$ and $E(\mathbf{Y};q)$.

We will use the semisimple group of adjoint type $G_2$ as an example. The command `G=coxgroup(:G,2)` selects this group. One can instead choose `coxgroup(:A,2)`, `coxgroup(:B,2)`, and so on. Alternatively, one can select these groups using `rootdatum(:pgl,3)` or `rootdatum(:so,5)`, or non-semisimple groups such as `rootdatum(:gl,2)`.

Then `EX(G,g,n)` returns $E(\mathbf{X};q)$, and `EY(G,g,n)` returns $E(\mathbf{Y};q)$.

For example:

```julia
julia> EX(G,0,3)
Pol{BigInt}: q⁸+6q⁷+20q⁶+58q⁵+180q⁴+58q³+20q²+6q+1

julia> EY(G,0,3)
Pol{BigInt}: q⁸+6q⁷+19q⁶+45q⁵+99q⁴
```

![](https://github.com/baileywhitbread/CharacterVarieties.jl/blob/main/animation.gif)




## To do
- Speed up $E$-polynomial calculations:
```julia
julia> @time EX(coxgroup(:F,4),1,1)
61.184813 seconds (227.91 M allocations: 18.563 GiB, 8.07% gc time, 88.52% compilation time)
```



- There is a large, unexplained leap in computation time between $`\mathrm{SO}_9`$ and $`\mathrm{SO}_{11}`$:
```julia
julia> @time EX(coxgroup(:B,3),0,3)
  0.622608 seconds (1.15 M allocations: 88.349 MiB, 96.49% compilation time)

julia> @time EX(coxgroup(:B,4),0,3)
  0.580247 seconds (7.36 M allocations: 715.632 MiB, 12.29% gc time, 28.14% compilation time)

julia> @time EX(coxgroup(:B,5),0,3)
127.007951 seconds (2.15 G allocations: 218.391 GiB, 17.85% gc time, 0.44% compilation time)
```

- People who understand Julia tell me I should reduce the allocation sizes...

- Add real-time calculations of $E$-polynomials to [baileywhitbread.com](https://www.baileywhitbread.com).




## Further directions

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
This is known in one case (because $H(\mathbf{X};q,t)$ is explicitly known): $G=\mathrm{GL}_2$, $n=1$ and $C_1$ is the conjugacy class of 
```math
\left(\begin{smallmatrix}-1 & \\ & -1 \end{smallmatrix}\right).
```
When $G=\mathrm{GL}_d$ and the $C_i$ are semisimple, an unproven formula for $H(\mathbf{X};q,t)$ was given in [HLRV](https://projecteuclid.org/journals/duke-mathematical-journal/volume-160/issue-2/Arithmetic-harmonic-analysis-on-character-and-quiver-varieties/10.1215/00127094-1444258.full). One could use the specialisations $E(\mathbf{X};q)$ and $E(\mathbf{Y};q)$ to try and understand the mixed Hodge polynomial $H(\mathbf{X};q,t)$.

The formulas for $E(\mathbf{X};q)$ and $E(\mathbf{Y};q)$ in [HLRV](https://projecteuclid.org/journals/duke-mathematical-journal/volume-160/issue-2/Arithmetic-harmonic-analysis-on-character-and-quiver-varieties/10.1215/00127094-1444258.full) are very different to ours. When $G=\mathrm{GL}_d$, the coefficients of $`\#\mathbf{Y}(\mathbb{F}_q)`$ are non-negative by the work of [HLRV](https://projecteuclid.org/journals/duke-mathematical-journal/volume-160/issue-2/Arithmetic-harmonic-analysis-on-character-and-quiver-varieties/10.1215/00127094-1444258.full). One could use `CharacterVarieties.jl` to search for $E(\mathbf{Y};q)$ with negative coefficients, or other interesting behavior.
