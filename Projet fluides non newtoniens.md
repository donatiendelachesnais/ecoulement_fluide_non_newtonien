## Equations

<span style="color: blue;">(1)</span>
$-\nabla.(\mu \nabla w) = \beta$

<span style="color: blue;">(2)</span> 
$\mu = \mu_T |\nabla w|^{n-1} = \mu_T (|\nabla w|^2)^{\frac{n-1}{2}}$

<span style="color: blue;">(3)</span>
$\mu_T(T) = \mu_0 \exp(-\alpha (T - T_0))$

<span style="color: blue;">(4)</span>
$-k \Delta T = \mu |\nabla w|^2$

## Formulation faible du problème (alpha = 0)

A partir du problème  <span style="color: blue;">(1)-(2)</span>, on obtient :
$$
\int_\Omega \mu(|\nabla w|^2)
\nabla w\cdot\nabla\varphi
=
\beta\int_\Omega\varphi
\qquad \forall\varphi\in H_D^1(\Omega).
$$


## Point fixe sans dépendance de la température
Calculer $\mu_{k+1} = \mu(|\nabla w_k|)$

et avec dépendance de la température :
$\mu_{k+1} = \mu_T(T_k) |\nabla w_k|^{n-1}$

puis résoudre :
<span style="color: blue;">(5)</span>
$-\nabla.(\mu_k \nabla w_{k+1}) = \beta$

La forme faible s’écrit :
$$
\int_{\Omega} \mu_{k+1} \, \nabla w_{k+1} \cdot \nabla v \, d\Omega
= \int_{\Omega} \beta v \, d\Omega,

\qquad T=100 \text{ sur } D, \quad T=0 \text{ sur } C_1,\quad \frac{\partial T}{\partial n}=0 \text{ sur } C_2.
$$

Avec la température, la forme faible s'écrit :

$$
\int_{\Omega} \nabla T \cdot \nabla v \, d\Omega
- \int_{\Omega} f v \, d\Omega
+ \int_{C_2} \mathrm{q}_{C_2} \, v \, ds = 0,

\qquad T=100 \text{ sur } D, \quad T=0 \text{ sur } C_1, \quad \frac{\partial T}{\partial n}=0 \text{ sur } C_2.
$$


## Méthode de Newton

Nous avons précédemment obtenu la formulation faible du problème initial, qui s'écrit :

$$
\int_\Omega
\mu(|\nabla w|^2)\nabla w\cdot\nabla\varphi
=
\beta\int_\Omega\varphi,
\qquad \forall\varphi\in H_D^1(\Omega).
$$

Cette formulation reste non linéaire puisque la viscosité $\mu$ dépend de la solution $w$, et plus précisément de son gradient $\nabla w$. Nous allons donc utiliser la méthode de Newton afin de linéariser le problème autour d'une approximation $w^k$ de la solution.

### Principe de la méthode de Newton

À l'itération $k$, on cherche une nouvelle approximation $w^{k+1}$ sous la forme :

$$
w^{k+1}=w^k+h,
$$

où $w^k$ est l'approximation actuelle, supposée connue, et $h$ est la correction que l'on cherche à déterminer.

Par conséquent, on a :

$$
\nabla w^{k+1}=\nabla w^k+\nabla h.
$$

Pour appliquer la méthode de Newton, nous allons linéariser le terme non linéaire de la formulation faible, à savoir :

$$
F(w)=\mu(|\nabla w|^2)\nabla w.
$$

L'idée consiste à effectuer un développement de Taylor au premier ordre de cette fonction autour de $w^k$, dans la direction de la correction $h$. On écrit ainsi :

$$
F(w^k+h)\approx F(w^k)+DF(w^k)[h],
$$

où $DF(w^k)[h]$ désigne la dérivée de $F$ au point $w^k$, dans la direction $h$.

Il nous faut donc calculer cette dérivée directionnelle.

### Calcul de la dérivée directionnelle

On considère la fonction :

$$
F(w)=\mu(|\nabla w|^2)\nabla w.
$$

Cette expression est le produit d'une fonction scalaire, $\mu(|\nabla w|^2)$, par le vecteur $\nabla w$. Pour calculer sa dérivée directionnelle, nous allons utiliser la règle de dérivation d'un produit.

On obtient :

$$
\begin{aligned}
DF(w^k)[h]
={}&D\left(\mu(|\nabla w|^2)\right)_{w^k}[h]
\,\nabla w^k\\
&+\mu(|\nabla w^k|^2)\nabla h.
\end{aligned}
$$

Il reste à calculer la dérivée de la fonction scalaire $\mu(|\nabla w|^2)$.

Pour cela, on rappelle que :

$$
|\nabla w|^2=\nabla w\cdot\nabla w.
$$

En dérivant cette expression dans la direction $h$, on obtient :

$$
\begin{aligned}
D(|\nabla w|^2)_{w^k}[h]
&=D(\nabla w\cdot\nabla w)_{w^k}[h]\\
&=\nabla h\cdot\nabla w^k
+\nabla w^k\cdot\nabla h\\
&=2\nabla w^k\cdot\nabla h.
\end{aligned}
$$

Puisque $\mu$ est une fonction de la variable scalaire $|\nabla w|^2$, on utilise la règle de dérivation en chaîne. Ainsi :

$$
\begin{aligned}
D\left(\mu(|\nabla w|^2)\right)_{w^k}[h]
&=\mu'(|\nabla w^k|^2)
D(|\nabla w|^2)_{w^k}[h]\\
&=2\mu'(|\nabla w^k|^2)
(\nabla w^k\cdot\nabla h).
\end{aligned}
$$

En remplaçant la dérivée directionnelle de $\mu$ calculée ci-dessus dans celle de $F$, on obtient :

$$
\boxed{
\begin{aligned}
DF(w^k)[h]
={}&\mu(|\nabla w^k|^2)\nabla h\\
&+2\mu'(|\nabla w^k|^2)
(\nabla w^k\cdot\nabla h)\nabla w^k.
\end{aligned}
}
$$

Nous avons donc calculé la dérivée directionnelle du terme non linéaire $F$, ce qui va nous permettre d'écrire le problème linéarisé de Newton.

### Loi régularisée utilisée dans le code

Dans `methode_newton.edp`, la loi (2) est régularisée afin d'éviter une singularité lorsque $\nabla w=0$. En posant $s=|\nabla w|^2$ et $q_\epsilon=\sqrt{s+\epsilon^2}$, le code utilise
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
L'exposant $n-3$ provient donc de la dérivation de la puissance $(s+\epsilon^2)^{(n-1)/2}$. Dans le code, `gradWNorm` vaut $q_\epsilon$, `muW` vaut $\mu_\epsilon(s)$ et `dMuDs` vaut $\mu_\epsilon'(s)$. Le Jacobien contient $2\,\texttt{dMuDs}$, ce qui donne le facteur $2\mu'(s)$ de la formule générale.

Les paramètres utilisés sont $\beta=200$, $\mu_T=50$, $n=1$ et $\epsilon=10^{-8}$. Avec $n=1$, $\mu_\epsilon(s)=\mu_T$ et $\mu_\epsilon'(s)=0$ : le terme contenant `dMuDs` s'annule et le problème en $w$ est linéaire. Pour $n\ne1$, la formule générale de Newton s'applique avec cette loi régularisée.

En notant $s_k=|\nabla w^k|^2$ et $q_k=\sqrt{s_k+\epsilon^2}$, le système de Newton assemblé par le code s'écrit explicitement
$$
\begin{aligned}
&\int_\Omega \mu_T q_k^{\,n-1}\nabla h\cdot\nabla\varphi\\
&+\int_\Omega \mu_T(n-1)q_k^{\,n-3}
(\nabla w^k\cdot\nabla h)(\nabla w^k\cdot\nabla\varphi)\\
&=\beta\int_\Omega\varphi
-\int_\Omega \mu_T q_k^{\,n-1}\nabla w^k\cdot\nabla\varphi.
\end{aligned}
$$

### Construction du problème linéarisé de Newton

On rappelle que la formulation faible initiale est :

$$
\int_\Omega F(w)\cdot\nabla\varphi
=
\beta\int_\Omega\varphi,
\qquad \forall\varphi\in H_D^1(\Omega).
$$

À l'itération $k$, on remplace $w$ par $w^{k+1}=w^k+h$, puis on utilise le développement de Taylor au premier ordre :

$$
F(w^{k+1})\approx F(w^k)+DF(w^k)[h].
$$

La formulation faible devient alors :

$$
\int_\Omega
\left(F(w^k)+DF(w^k)[h]\right)
\cdot\nabla\varphi
=
\beta\int_\Omega\varphi.
$$

En utilisant la linéarité de l'intégrale, on sépare les deux termes :

$$
\begin{aligned}
&\int_\Omega F(w^k)\cdot\nabla\varphi\\
&+\int_\Omega DF(w^k)[h]\cdot\nabla\varphi
=\beta\int_\Omega\varphi.
\end{aligned}
$$

On cherche maintenant à isoler les termes contenant la correction $h$. Pour cela, on soustrait le terme $\displaystyle\int_\Omega F(w^k)\cdot\nabla\varphi$ des deux côtés de l'égalité. On obtient :

$$
\begin{aligned}
\int_\Omega DF(w^k)[h]\cdot\nabla\varphi
={}&\beta\int_\Omega\varphi\\
&-\int_\Omega F(w^k)\cdot\nabla\varphi.
\end{aligned}
$$

Nous remplaçons maintenant $DF(w^k)[h]$ par son expression calculée précédemment, ainsi que $F(w^k)$ par sa définition.

On obtient finalement :

$$
\boxed{
\begin{aligned}
&\int_\Omega
\mu(|\nabla w^k|^2)
\nabla h\cdot\nabla\varphi\\
&+2\int_\Omega
\mu'(|\nabla w^k|^2)
(\nabla w^k\cdot\nabla h)
(\nabla w^k\cdot\nabla\varphi)\\
&=\beta\int_\Omega\varphi
-\int_\Omega
\mu(|\nabla w^k|^2)
\nabla w^k\cdot\nabla\varphi,
\qquad \forall\varphi\in H_D^1(\Omega).
\end{aligned}
}
$$

On retrouve ainsi la formulation variationnelle (6) donnée dans l'énoncé.

À chaque itération, le problème consiste donc à déterminer la correction $h$ à partir d'un problème linéaire, puisque les coefficients $\mu(|\nabla w^k|^2)$ et $\mu'(|\nabla w^k|^2)$ sont calculés à partir de l'approximation connue $w^k$.

Une fois $h$ obtenue, on met à jour la solution par :

$$
\boxed{w^{k+1}=w^k+h.}
$$

On répète ensuite cette procédure jusqu'à ce que la correction soit suffisamment petite, ce qui indique que l'algorithme a convergé vers une solution du problème non linéaire.


### Lien entre les formulations faible et forte

Posons
$$
F_k=\mu(|\nabla w^k|^2)\nabla w^k,
\qquad
M_k=M(\nabla w^k)
=\mu(|\nabla w^k|^2)I
+2\mu'(|\nabla w^k|^2)\,\nabla w^k\otimes\nabla w^k.
$$
Avec la convention $(a\otimes a)b=a(a\cdot b)$, le calcul de la dérivée
directionnelle donne bien $DF(w^k)[h]=M_k\nabla h$. La formulation (6) devient
alors
$$
\int_\Omega M_k\nabla h\cdot\nabla\varphi
=\beta\int_\Omega\varphi-\int_\Omega F_k\cdot\nabla\varphi.
$$

En intégrant par parties, le terme $-\int_\Omega F_k\cdot\nabla\varphi$
contribue $+\int_\Omega(\nabla\cdot F_k)\varphi$ (à un terme de bord près).
On obtient donc, dans $\Omega$,
$$
-\nabla\cdot\left(M(\nabla w^k)\nabla h\right)
=\nabla\cdot\left(\mu(|\nabla w^k|^2)\nabla w^k\right)+\beta.
$$
Le signe indiqué dans l'équation proposée est donc erroné : il faut un signe
$+$ devant la divergence de $\mu(|\nabla w^k|^2)\nabla w^k$. Sur la partie
de Neumann, la condition naturelle correspondante porte sur le flux total,
$(M_k\nabla h+F_k)\cdot n=0$; les termes de bord s'annulent aussi pour les
fonctions tests nulles sur la partie de Dirichlet.



## Problème couplé

On calcule les suites $w_l$ et $T_l$.

## Paramètres et autres formules de `methode_newton.edp`

Le domaine numérique est délimité par les arcs du cercle unité
$$
C_1:\ (x,y)=(\cos t,\sin t),\quad 0\le t\le\pi,
\qquad
C_2:\ (x,y)=(\cos t,\sin t),\quad \pi\le t\le2\pi,
$$
et par le cercle intérieur
$$
D:\ (x,y)=(0.3+0.5\cos t,\ 0.5\sin t),\quad 0\le t\le2\pi.
$$
Le maillage est construit par `buildmesh(C1(25)+C2(25)+D(-30))` avec des éléments finis $P_1$. La correction vérifie $h=0$ sur $C_1$, $C_2$ et $D$; l'itération initiale est $w^0=0$.

Le critère d'arrêt du code est la norme $L^2$ de la correction :
$$
\|h\|_{L^2(\Omega)}
=\left(\int_\Omega h^2\right)^{1/2}<10^{-6},
$$
avec au plus $50$ itérations de Newton. Les paramètres `tolerance` et `maxIterations` déclarés pour la méthode du point fixe ne sont pas utilisés par cette boucle.

Pour la température, le code utilise $k=1$ et définit
$$
f=\frac{\mu}{k}|\nabla w|^2
=\frac{\mu}{k}\left((\partial_x w)^2+(\partial_y w)^2\right).
$$
Il résout la formulation faible déjà donnée plus haut, avec $T=0$ sur $C_1$, $T=100$ sur $D$ et `fluxC2` nul sur $C_2$. Le champ `mu` est initialisé à $\mu_T$ et n'est pas recalculé après la résolution de $w$; cela correspond à la loi utilisée ici puisque $n=1$. Si $n$ est modifié, il faudra recalculer `mu` à partir de $w$ pour représenter la viscosité non linéaire dans $f$.

Enfin, si `fluxC2` désigne le flux normal sortant $\partial T/\partial n$, la formulation faible usuelle comporte $-\int_{C_2}\mathrm{fluxC2}\,v\,ds$. Le code assemble le terme avec un signe $+$; cela n'affecte pas le calcul présent car `fluxC2=0`.