# epolys
We have computed E-polynomials of character varieties associated to general reductive groups. 

Using Jean Michel's package 'Chevie', these polynomials are computable in Julia. 

I have written a Julia script that computes these E-polynomials. 

This works on Julia v1.10.2 using the Chevie version dated April 12th 2024. 

The setting is as follows: 

Let $G$ be a connected split reductive group over $\mathbb{F}_q$ with connected centre $Z$ and split maximal torus $T$. Fix integers $g\geq 0$ and $n\geq 1$, and select a strongly regular elements $S_1,\ldots,S_n$ in $T$ that are 'generic' in the sense of our paper. Let $C_1,\ldots,C_n$ be their conjugacy classes. Then the character variety is the GIT quotient
```math
\mathbf{X} = \bigg\{(A_1,B_1,\ldots,A_g,B_g,Y_1,\ldots,Y_n)\in G^{2g}\times C_1\times \cdots\times C_n\ \bigg|\ [A_1,B_1]\cdots[A_g,B_g]Y_1\cdots Y_n = 1\bigg\}\bigg/\!\!\!\!\bigg/G
```
where the action is simultaneous conjugation. 

The Julia script computes $|\mathbf{X}(\mathbb{F}_q)|$. 
