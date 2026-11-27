    ## Equations

Equation 1 : 
$-\nabla.(\mu \nabla w) = \beta$

Eq. 2 : 
$\mu = \mu_T |\nabla w|^{n-1} = \mu_T (|\nabla w|^2)^{(n-1)/2}$

Eq. 3 :
$\mu_T(T) = \mu_0 \exp(-\alpha (T - T_0))$

Eq. 4 :
$-k \Delta T = \mu |\nabla w|^2$

## Point fixe sans dépendance de la température (alpha = 0)
Calculer $\mu_{k+1} = \mu(|\nabla w_k|)$

et avec dépendance de la température :
$\mu_{k+1} = \mu_T(T_k) |\nabla w_k|^{n-1}$

puis résoudre :
$-\nabla.(\mu_k \nabla w_{k+1}) = \beta$

## Méthode de Newton sans dépendance de la température (alpha = 0)
On note $w_k$ l'approximation de $w$ à l'itération $k$.

Montrons que si $h$ est la correction $w_{k+1} = w_k + h$, alors $h$ est solution du problème variationnel suivant :

Eq. 6:
$$\int_\Omega \left( \mu(|\nabla w^k|^2)\nabla h \cdot \nabla \varphi + 2\mu'(|\nabla w^k|^2)(\nabla w^k \cdot \nabla \varphi)(\nabla w^k \cdot \nabla h) \right) = \beta \int_\Omega \varphi - \int_\Omega \mu(|\nabla w^k|^2)\nabla w^k \cdot \nabla \varphi, \quad \forall \varphi \in H_D^1(\Omega)$$

Posons $s(w)=|\nabla w|^2$ et définissons le résidu faible
$$
F(w)(\varphi)
= \int_\Omega \mu(s(w))\,\nabla w\cdot\nabla\varphi
- \beta\int_\Omega \varphi.
$$
La formulation faible de l'équation $-\nabla\cdot(\mu(|\nabla w|^2)\nabla w)=\beta$ consiste à chercher $w$ tel que $F(w)(\varphi)=0$ pour tout $\varphi\in H_D^1(\Omega)$.

À l'itération $k$, on cherche une correction $h$ telle que $w_{k+1}=w_k+h$. La méthode de Newton remplace le résidu par son approximation linéaire en $w_k$ :
$$
F(w_k+h)(\varphi)\simeq F(w_k)(\varphi)+DF(w_k)[h](\varphi)=0,
$$
où $DF(w_k)[h]$ désigne la dérivée de $F$ en $w_k$ dans la direction $h$.

Pour calculer cette dérivée, posons $s_k=|\nabla w_k|^2$ et considérons $w_k+t h$. Alors
$$
\left.\frac{d}{dt}s(w_k+t h)\right|_{t=0}
=2\nabla w_k\cdot\nabla h.
$$
Dans le code, la loi de viscosité est régularisée pour éviter une singularité lorsque $\nabla w=0$. En posant $s=|\nabla w|^2$ et $q_\epsilon=\sqrt{s+\epsilon^2}$, elle devient
$$
\mu_\epsilon(s)=\mu_T q_\epsilon^{\,n-1}
=\mu_T(s+\epsilon^2)^{(n-1)/2}.
$$
Sa dérivée par rapport à $s$ est
$$
\begin{aligned}
\mu_\epsilon'(s)
&=\frac{d}{ds}\left[\mu_T(s+\epsilon^2)^{(n-1)/2}\right]\\
&=\frac{\mu_T(n-1)}{2}(s+\epsilon^2)^{(n-3)/2}\\
&=\frac{\mu_T(n-1)}{2}q_\epsilon^{\,n-3}.
\end{aligned}
$$
L'exposant $n-3$ vient donc de la dérivation de la puissance $(s+\epsilon^2)^{(n-1)/2}$. Dans le code, `gradWNorm` représente $q_\epsilon$, `muW` représente $\mu_\epsilon(s)$ et `dMuDs` représente $\mu_\epsilon'(s)$. Le terme du Jacobien utilise $2\,\texttt{dMuDs}$, conformément au facteur $2\mu'(s)$ de la dérivée directionnelle.

Les paramètres actuellement utilisés sont $\beta=200$, $\mu_T=50$, $n=1$ et $\epsilon=10^{-8}$. Pour $n=1$, on a $\mu_\epsilon(s)=\mu_T$ et $\mu_\epsilon'(s)=0$ : le terme contenant `dMuDs` s'annule et le problème en $w$ est linéaire. Pour $n\ne1$, la formule générale de Newton ci-dessus s'applique avec cette loi régularisée.

En notant $s_k=|\nabla w_k|^2$ et $q_k=\sqrt{s_k+\epsilon^2}$, le système assemblé par le code s'écrit explicitement
$$
\begin{aligned}
&\int_\Omega \mu_T q_k^{\,n-1}\nabla h\cdot\nabla\varphi\\
&+\int_\Omega \mu_T(n-1)q_k^{\,n-3}
(\nabla w_k\cdot\nabla h)(\nabla w_k\cdot\nabla\varphi)\\
&=\beta\int_\Omega\varphi
-\int_\Omega \mu_T q_k^{\,n-1}\nabla w_k\cdot\nabla\varphi.
\end{aligned}
$$

Par la règle du produit et la règle de la chaîne, (ici, $\mu'$ désigne la dérivée de $\mu$ par rapport à son argument scalaire $s$)
$$
\left.\frac{d}{dt}\left(\mu(s(w_k+t h))\nabla(w_k+t h)\right)\right|_{t=0}
=2\mu'(s_k)(\nabla w_k\cdot\nabla h)\nabla w_k
+\mu(s_k)\nabla h.
$$
En prenant le produit scalaire avec $\nabla\varphi$ et en intégrant sur $\Omega$, on obtient
$$
DF(w_k)[h](\varphi)
=\int_\Omega\left(
\mu(s_k)\nabla h\cdot\nabla\varphi
+2\mu'(s_k)(\nabla w_k\cdot\nabla h)(\nabla w_k\cdot\nabla\varphi)
\right).
$$
Comme $F(w_k)(\varphi)=\int_\Omega\mu(s_k)\nabla w_k\cdot\nabla\varphi-\beta\int_\Omega\varphi$, la condition $DF(w_k)[h](\varphi)=-F(w_k)(\varphi)$ donne bien l'équation annoncée. Enfin, si $w$ est imposée sur le bord, la correction $h$ a une trace nulle sur cette partie du bord, d'où $h\in H_D^1(\Omega)$.


#### Formulation faible

Posons $A_k=\mu(|\nabla w_k|^2)\nabla w_k$. Avec
$$
M(\boldsymbol{\xi})
=\mu(|\boldsymbol{\xi}|^2)I
+2\mu'(|\boldsymbol{\xi}|^2)\,\boldsymbol{\xi}\otimes\boldsymbol{\xi},
$$
on a
$$
M(\nabla w_k)\nabla h\cdot\nabla\varphi
=\mu(|\nabla w_k|^2)\nabla h\cdot\nabla\varphi
+2\mu'(|\nabla w_k|^2)(\nabla w_k\cdot\nabla h)(\nabla w_k\cdot\nabla\varphi).
$$
L'équation (6) devient donc
$$
\int_\Omega M(\nabla w_k)\nabla h\cdot\nabla\varphi
=\beta\int_\Omega\varphi-\int_\Omega \mu(|\nabla w_k|^2)\nabla w_k\cdot\nabla\varphi.
$$
En supposant que $\varphi$ s'annule sur le bord (ou que les termes de bord s'annulent par les conditions aux limites), l'intégration par parties donne
$$
\int_\Omega M(\nabla w_k)\nabla h\cdot\nabla\varphi
=\int_\Omega\left[-\nabla\cdot\left(M(\nabla w_k)\nabla h\right)\right]\varphi,
\qquad
-\int_\Omega A_k\cdot\nabla\varphi
=\int_\Omega\left[\nabla\cdot A_k\right]\varphi.
$$
Par conséquent, la forme forte équivalente à (6) est
$$
-\nabla\cdot\left(M(\nabla w_k)\nabla h\right)
=\nabla\cdot\left(\mu(|\nabla w_k|^2)\nabla w_k\right)+\beta,
\qquad\text{dans }\Omega.
$$

L'équation de l'image indique cependant
$$
-\nabla\cdot\left(M(\nabla w^k)\nabla h\right)
=-\nabla\cdot\left(\mu\nabla w^k\right)+\beta,
\qquad \text{dans }\Omega,
$$
avec
$$
M(\boldsymbol{\xi})
=\mu(|\boldsymbol{\xi}|^2)I
+2\mu'(|\boldsymbol{\xi}|^2)\,\boldsymbol{\xi}\otimes\boldsymbol{\xi}.
$$
Le désaccord apparaît au second terme : dans (6), on a $-\int_\Omega A_k\cdot\nabla\varphi$, qui devient $+\int_\Omega(\nabla\cdot A_k)\varphi$ par intégration par parties. Il faut donc $+\nabla\cdot A_k$ dans la forme forte, et non $-\nabla\cdot A_k$. Pour obtenir le signe de l'image, le terme correspondant dans (6) devrait être $+\int_\Omega A_k\cdot\nabla\varphi$.

#### De la forme forte à la forme faible

Réciproquement, partons de la forme forte cohérente avec (6) :
$$
-\nabla\cdot\left(M(\nabla w_k)\nabla h\right)
=\nabla\cdot A_k+\beta
\qquad\text{dans }\Omega.
$$
Multiplions cette équation par une fonction test $\varphi$ et intégrons sur $\Omega$. La formule de Green donne
$$
\int_\Omega M(\nabla w_k)\nabla h\cdot\nabla\varphi
-\int_{\partial\Omega}\varphi\,\left(M(\nabla w_k)\nabla h\right)\cdot n\,dS
=\beta\int_\Omega\varphi
-\int_\Omega A_k\cdot\nabla\varphi
+\int_{\partial\Omega}\varphi\,A_k\cdot n\,dS.
$$
Si $\varphi$ s'annule sur le bord, les deux intégrales de bord sont nulles. On obtient alors
$$
\int_\Omega M(\nabla w_k)\nabla h\cdot\nabla\varphi
=\beta\int_\Omega\varphi-\int_\Omega A_k\cdot\nabla\varphi,
$$
qui est précisément l'équation (6), puisque le développement de $M(\nabla w_k)\nabla h\cdot\nabla\varphi$ donne les deux termes du membre de gauche de cette équation. Si $\varphi$ ne s'annule pas sur toute la frontière, les termes de bord doivent être conservés ou annulés par les conditions aux limites.

## Formules d'implémentation complémentaires

### Domaine et conditions aux limites

Le bord extérieur est le cercle unité, découpé en deux arcs,
$$
C_1:\ (x,y)=(\cos t,\sin t),\quad 0\le t\le\pi,
\qquad
C_2:\ (x,y)=(\cos t,\sin t),\quad \pi\le t\le2\pi,
$$
et le bord intérieur est le cercle de rayon $0.5$ centré en $(0.3,0)$,
$$
D:\ (x,y)=(0.3+0.5\cos t,\ 0.5\sin t),\quad 0\le t\le2\pi.
$$
Le maillage est construit avec l'appel `buildmesh(C1(25)+C2(25)+D(-30))` et utilise des éléments finis $P_1$. La correction de Newton satisfait $h=0$ sur $C_1$, $C_2$ et $D$. Avec l'initialisation $w^0=0$, ces conditions maintiennent aussi $w=0$ sur ces bords.

### Critère d'arrêt de Newton

Le code calcule à chaque itération la norme $L^2$ de la correction,
$$
\|h\|_{L^2(\Omega)}
=\left(\int_\Omega h^2\right)^{1/2},
$$
et arrête les itérations si elle est inférieure à $10^{-6}$, ou après au plus $50$ itérations. Les variables `tolerance` et `maxIterations` déclarées plus haut concernent l'ancienne méthode du point fixe et ne contrôlent pas cette boucle de Newton.

### Équation de température

L'équation $-k\Delta T=\mu|\nabla w|^2$ est divisée par $k$ dans le code, qui définit
$$
f=\frac{\mu}{k}|\nabla w|^2
=\frac{\mu}{k}\left((\partial_x w)^2+(\partial_y w)^2\right).
$$
La forme faible assemblée par le code est
$$
\int_\Omega\nabla T\cdot\nabla v
-\int_\Omega f v
+\int_{C_2}q_{C_2}v\,ds=0,
\qquad
T=0\text{ sur }C_1,\quad T=100\text{ sur }D,\quad q_{C_2}=0.
$$
Dans ce fichier, $k=1$ et le champ `mu` vaut $\mu_T$; avec $n=1$, cela correspond à la loi de viscosité utilisée pour $w$. Si $n$ est modifié, le champ `mu` employé dans $f$ doit aussi être recalculé à partir de $w$ pour représenter la viscosité non linéaire.

Remarque de signe : si $q_{C_2}$ désigne le flux normal sortant $\partial T/\partial n$, la formulation faible usuelle comporte $-\int_{C_2}q_{C_2}v\,ds$. Le code utilise le signe $+$; ici cela n'affecte pas la solution puisque $q_{C_2}=0$.