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

En remplaçant cette expression dans le calcul de $DF(w^k)[h]$, on obtient :

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

## Problème couplé

On calcule les suites $w_l$ et $T_l$.
