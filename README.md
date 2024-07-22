# CharacterVarieties.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://baileywhitbread.github.io/CharacterVarieties.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://baileywhitbread.github.io/CharacterVarieties.jl/dev/)
[![Build Status](https://github.com/baileywhitbread/CharacterVarieties.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/baileywhitbread/CharacterVarieties.jl/actions/workflows/CI.yml?query=branch%3Amaster)


This package computes E-polynomials of character varieties associated to general reductive groups in Julia using Jean Michel's package [Chevie](https://github.com/jmichel7/Chevie.jl). 

This works on Julia v1.10.2 using the Chevie version dated April 12th 2024. 

## Multiplicative character varieties
Let $G$ be a connected split reductive group over $\mathbb{F}_q$ with connected centre $Z$ and split maximal torus $T$. Fix integers $g\geq 0$ and $n\geq 1$, and select a strongly regular elements $S_1,\ldots,S_n$ in $T$ that are 'generic.' Let $C_1,\ldots,C_n$ be their conjugacy classes. 

The multiplicative character variety is the GIT quotient
```math
\mathbf{X} = \bigg\{(A_1,B_1,\ldots,A_g,B_g,Y_1,\ldots,Y_n)\in G^{2g}\times C_1\times \cdots\times C_n\ \bigg|\ [A_1,B_1]\cdots[A_g,B_g]Y_1\cdots Y_n = 1\bigg\}\bigg/\!\!\!\!\!\bigg/G
```
where the action is simultaneous conjugation. 

This package computes the E-polynomial $E(\mathbf{X};q)$ via our formula for $\#\mathbf{X}(\mathbb{F}_q)$.  

## Additive character varieties
Let $\mathfrak{g}$ be the Lie algebra of $G$ and let $\mathfrak{t}$ be the Lie algebra of $T$. Select regular semisimple elements $s_1,\ldots,s_n$ in $\mathfrak{t}$ that are 'generic.' Let $O_1,\ldots,O_n$ be their adjoint orbits. 

The additive character variety is the GIT quotient
```math
\mathbf{Y} = \bigg\{(A_1,B_1,\ldots,A_g,B_g,Y_1,\ldots,Y_n)\in \mathfrak{g}^{2g}\times O_1\times \cdots\times O_n\ \bigg|\ [A_1,B_1]+\cdots+[A_g,B_g] + Y_1+ \cdots + Y_n = 0\bigg\}\bigg/\!\!\!\!\!\bigg/G
```
where the action is simultaneous conjugation (i.e., the adjoint action). 

This package computes the E-polynomial $E(\mathbf{Y};q)$ via our formula for $`\#\mathbf{Y}(\mathbb{F}_q)`$.  

## Mixed Hodge polynomials
Associated to $\mathbf{X}$ is the (compactly supported) mixed Hodge polynomial
```math
H(\mathbf{X};x,y,t) = \sum_{i,j,k} h^{i,j,k} x^i y^j t^k.
```
Proving that $`\#\mathbf{X}(\mathbb{F}_q)`$ is a polynomial in $q$ implies $H(\mathbf{X};x,y,t)$ depends only on the product $xy=:q$ and $t$, and that
```math
H(\mathbf{X};q,-1) = E(\mathbf{X};q).
```
That is, computing $`\#\mathbf{X}(\mathbb{F}_q)`$ gives us a specialisation of $H(\mathbf{X};q,t)$. 

## Pure parts
Another important specialisation of $H(\mathbf{X};q,t)$ is obtained in the following manner. Contained in the cohomology ring $H^\ast(\mathbf{X})$ is an important subring called the pure subring $H^\ast_\mathrm{pure}(\mathbf{X})$. Denote by $PH(\mathbf{X};u)$ the Poincar\'{e} polynomial of the pure subring. Then there is the specialisation $H(\mathbf{X})\mapsto PH(\mathbf{X})$ by setting all terms to zero except the monomials in $qt^2$ (for instance, if $H(\mathbf{X};q,t)=qt^2 + q^2t^4 + q + qt + 1$ then $PH(\mathbf{X};q)=q+q^2$). We say a variety's cohomology is 'pure' if $H=PH$. 


## Non-negative coefficients
When $G=\mathrm{GL}_n$, it was proven in [HLRV](https://projecteuclid.org/journals/duke-mathematical-journal/volume-160/issue-2/Arithmetic-harmonic-analysis-on-character-and-quiver-varieties/10.1215/00127094-1444258.full) that the coefficients of $`\#\mathbf{Y}(\mathbb{F}_q)`$ are non-negative. This was proven by relating the additive character variety to a quiver variety (which are known to have pure cohomology due to their symplectic geometry). Moreover, it was conjectured the polynomials $\mathbf{X}$ and $\mathbf{X}$ are closely related, in the sense that

```math
PH(\mathbf{X};q) = E(\mathbf{Y};q).
```

The conjecture has been proven in one narrow case: $G=\mathrm{GL}_2$, $n=1$ and $C_1$ contains $`\left(\begin{smallmatrix}-1 & \\ & -1 \end{smallmatrix}\right)`$. An unproven conjectural formula is given for $G=\mathrm{GL}_n$ in [HLRV](https://projecteuclid.org/journals/duke-mathematical-journal/volume-160/issue-2/Arithmetic-harmonic-analysis-on-character-and-quiver-varieties/10.1215/00127094-1444258.full). 

## The idea
Use Julia to generate the specialisations $E(\mathbf{X};q)$ and $E(\mathbf{Y};q)$ and try to understand the mixed Hodge polynomial $H(\mathbf{X};q,t)$. 
