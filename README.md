# CharacterVarieties.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://baileywhitbread.github.io/CharacterVarieties.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://baileywhitbread.github.io/CharacterVarieties.jl/dev/)
[![Build Status](https://github.com/baileywhitbread/CharacterVarieties.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/baileywhitbread/CharacterVarieties.jl/actions/workflows/CI.yml?query=branch%3Amaster)


This package computes counting polynomials of character varieties associated to general reductive groups in Julia using Jean Michel's package [Chevie](https://github.com/jmichel7/Chevie.jl). 

This works on Julia v1.10.2 using the Chevie version dated April 12th 2024. 

## Multiplicative character varieties:
Let $G$ be a connected split reductive group over $\mathbb{F}_q$ with connected centre $Z$ and split maximal torus $T$. Fix integers $g\geq 0$ and $n\geq 1$, and select a strongly regular elements $S_1,\ldots,S_n$ in $T$ that are 'generic.' Let $C_1,\ldots,C_n$ be their conjugacy classes. 

The multiplicative character variety is the GIT quotient
```math
\mathbf{X} = \bigg\{(A_1,B_1,\ldots,A_g,B_g,Y_1,\ldots,Y_n)\in G^{2g}\times C_1\times \cdots\times C_n\ \bigg|\ [A_1,B_1]\cdots[A_g,B_g]Y_1\cdots Y_n = 1\bigg\}\bigg/\!\!\!\!\!\bigg/G
```
where the action is simultaneous conjugation. 

One can compute counting polynomials $E(\mathbf{X};q)$ via formulas for $|\mathbf{X}(\mathbb{F}_q)|$.  

## Additive character varieties:
Let $\mathfrak{g}$ be the Lie algebra of $G$ and let $\mathfrak{t}$ be the Lie algebra of $T$. Select regular semisimple elements $s_1,\ldots,s_n$ in $\mathfrak{t}$ that are 'generic.' Let $O_1,\ldots,O_n$ be their adjoint orbits. 

The additive character variety is the GIT quotient
```math
\mathbf{Y} = \bigg\{(A_1,B_1,\ldots,A_g,B_g,Y_1,\ldots,Y_n)\in \mathfrak{g}^{2g}\times O_1\times \cdots\times O_n\ \bigg|\ [A_1,B_1]+\cdots+[A_g,B_g] + Y_1+ \cdots + Y_n = 0\bigg\}\bigg/\!\!\!\!\!\bigg/G
```
where the action is simultaneous conjugation (i.e., the adjoint action). 

One can compute counting polynomials $E(\mathbf{Y};q)$ via formulas for $|\mathbf{Y}(\mathbb{F}_q)|$.  
