# CharacterVarieties.jl
We have computed E-polynomials of character varieties associated to general reductive groups. 

Using Jean Michel's package [Chevie](https://github.com/jmichel7/Chevie.jl), these polynomials are computable in Julia. 

This Julia script computes these E-polynomials. 

This works on Julia v1.10.2 using the Chevie version dated April 12th 2024. 

The setting is as follows: 

Let $G$ be a connected split reductive group over $\mathbb{F}_q$ with connected centre $Z$ and split maximal torus $T$. Fix integers $g\geq 0$ and $n\geq 1$, and select a strongly regular elements $S_1,\ldots,S_n$ in $T$ that are 'generic.' Let $C_1,\ldots,C_n$ be their conjugacy classes. 

The character variety is the GIT quotient
```math
\mathbf{X} = \bigg\{(A_1,B_1,\ldots,A_g,B_g,Y_1,\ldots,Y_n)\in G^{2g}\times C_1\times \cdots\times C_n\ \bigg|\ [A_1,B_1]\cdots[A_g,B_g]Y_1\cdots Y_n = 1\bigg\}\bigg/\!\!\!\!\bigg/G
```
where the action is simultaneous conjugation. 

The Julia script computes the E-polynomial $E(\mathbf{X};q) = |\mathbf{X}(\mathbb{F}_q)|$. 

Our formula for the E-polynomial is as follows:

Consider the collection of pairs $(L,\rho)$ where $L$ is an endoscopy group of $G$ and $\rho$ is an irreducible representation of $W(L)$, the Weyl group of $L$. Then $W$ acts on this collection of pairs. We call an orbit $[(L,\rho)]=[L,\rho]$ a type of $G$. 

For each type, define:
1. $m_\tau(q) := q^{|\Phi(G)^+|-|\Phi(L)^+|} \frac{|L(\mathbb{F}_q)|}{\tilde{\rho}(1)}$, where
   1. $\Phi(L)^+$ is the positive roots of $L$, and
   2. $\tilde{\rho}$ is the (principal) unipotent character of $L$ corresponding to the irreducible representation $\rho$ of $W(L)$.
3. $\gamma_\tau = \frac{\dim(\rho)^n |W|^{n-1} |[L]|  \nu(L)}{|W(L)|^{n-1}}$, where
   1. $|[L]|$ is the orbit of $L$ under the action of $W$ on the collection of endoscopy groups,
   2. $\nu(L) := \sum_{L'\supseteq L} \pi_0^{L'} \mu(L,L')$, where
      1. The sum is over all isolated endoscopy groups containing $L$,
      2. $\pi_0^{L'}$ is the number of connected components of the centre of $L'$, and
      3. $\mu$ is the Mobius function on the poset of endoscopy groups ordered by inclusion.

Then
```math
E(\mathbf{X};q) = \frac{|Z(\mathbb{F}_q)|^2}{|T(\mathbb{F}_q)|^n} \sum_{\tau=[L,\rho]} \gamma_\tau m_\tau(q)^{2g-2+n}.
```
