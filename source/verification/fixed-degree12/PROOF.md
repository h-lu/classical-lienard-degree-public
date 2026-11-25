# Finite certificate for the half-return comparisons

Let \(G_s(x)=F(sx)\) on \(x\geq0\) for \(s=\pm1\). The forward half-fields are \((\dot x,\dot y)=(y-G_s(x),-x)\). The change \((x,y,\tau)\mapsto(-x,y,-\tau)\) relates the negative side to the backward left half of the original field. Write \(w_s(h)\) for the first lower-axis hit magnitude from \((0,h)\) and \(D(h)=w_+(h)-w_-(h)\). Equal magnitudes give a periodic orbit. Returning arcs are ordered as graphs over \(y\), since \(\dot y=-x<0\) in the open half-plane. The half-returns increase strictly on intervals bounded by returning arcs. A transverse upper trap also gives returns at every smaller positive height when \(|F'(0)|<2\). A trapped nonreturning orbit must tend to the origin, where the polar angular speed is bounded above by \(-1+|F'(0)|/2+o(1)\). This excludes convergence to the origin within the open half-plane. Each periodic orbit meets the upper ray exactly once by transversality and order. Section 2 of the main article proves these return properties. Appendix SB of the electronic supplement gives the proof for the specified rational field and the fixed-scale matching estimates.

`rational.py` rounds interval bounds outwards to multiples of \(2^{-4096}\) using integer floor and ceiling. The mesh applies to the algebraic estimates for logarithms. Bernstein coefficients certify polynomial signs on closed rational intervals. The scale, profile and shooting constants are in `field.py`, with \(\eta=\lambda^2\) and \(E'(z)=z^5(1-z^2)^2(1-z^2/d^2)\). The fixed coefficients are in `inputs/fixed-diagnostic-field.json` and `outer.py` compares them with the independently evaluated rational formulas. Approximate numbers in the transcript are display roundings.

## 1. The backward outer passage

Write \(P(z)=F(sLz)/S\), \(u=z^2/2\), \(p(u)=P(\sqrt{2u})\), \(t=y/S\), and \(\epsilon=L^2/S^2\). In backward physical time the height increases while \(x>0\), and the exact scalar equation is

\[
 \epsilon u_t=p(u)-t.
\]

`hypotheses` in `outer.py` verifies the following bounds with \(K=10^7\), \(\rho=\lambda^4\), \(c=10^{-4}\), and \(\nu=\rho^2\).

- The unique global maximum \(M=p(u_m)\in(0.012,0.014)\) on \(0\leq u\leq2\) lies in \(1.18<\sqrt{2u_m}<1.22\). The derivative signs at the ends, \(p''<-0.1\) in this window, and \(p'<-0.01\) on \(1.22\leq z\leq2\) give uniqueness there. For \(z<1.18\), the leading even profile \(E(z)\) is increasing and \(P(d)>E(1.18)+\|P-E\|_{[0,2]}\). The same comparison at \(z=1.11\) gives \(P(z)<M-2c\) for \(0\leq z\leq1.11\).
- \(P>10^{-3}\) on \(1/2\leq z\leq1.22\), \(P>-10^{-3}\) on \(0\leq z\leq1.22\), \(P(1.22)>0.01\), and \(P(2)<-2\). The first two bounds use the increase of \(E\) up to \(d\) and its decrease from \(d\) to 1.22. The tail inequality \(P(z)<-10^{-3}z^{12}\) for \(z\geq2\) follows from the negative coefficients after translating \(P(z)+10^{-3}z^{12}\) by 2.
- On \(0\leq u\leq2\) one has \(|p|+2<K\). On \(z\geq1/2\) within that range one has \(|p''|<K\). The divided difference of \(p\) between two points of the full range has magnitude bounded by \(K\) whenever at least one of their \(z\) coordinates is at least \(1/2\). The conversion factor from the spatial divided difference is \(2/(z_1+z_2)\leq4\).

Let \(U(t)\) be the falling root of \(p(U)=t\) for \(-1\leq t\leq M-\rho\). On \([U(t)-\nu,2]\) one has \(-p'\geq\rho\). The bound on \(|p''|\) would give \(M-t\leq200K\rho^2<\rho\) if \(U(t)-u_m\leq20\rho\). Thus if \(U(t)-\nu\) is still in the maximum window, the curvature bound gives \(-p'(U(t)-\nu)>0.1(20\rho-\nu)>\rho\). To the right of that window use \(-p'>0.01\). In particular \(|U'|\leq1/\rho\).

The tube \(|u-U(t)|\leq\nu\) is invariant in increasing \(t\). Its boundary restoring speed is at least \(\rho\nu/\epsilon\), and \(\epsilon<\rho^4/8\). Consider the reference launch \(u=0,t=-1\), the direct launch \(u=0,t=-1/2\), or a section launch \(u=1/8,t=r=y/S\) with \(|y|\leq25\). Each enters the tube within a height increase of \(\rho\). The wall \(u=2\) points inward. Below \(U(t)-\nu\) on the falling branch the speed toward the tube is at least \(\rho\nu/\epsilon\) by integrating the slope over the \(\nu\) interval and using monotonicity. Below the falling branch during the launch layer the same speed bound follows from the bounds on \(P\) through \(z=1.22\), since \(U(t)-\nu\) lies to the right of \(z=1.22\) in that layer. The positive-height section launches stay beyond \(z=1/2\) because \(P(1/2)>10^{-3}\) and \(|r|+\rho<2\rho<10^{-3}\). Above the tube the slope bound gives speed at least \(\rho\nu/\epsilon\) toward it. Even accounting for root motion, traversing distance 2 requires less than \(\rho\) in height. Axis launches at negative \(t\) enter \(x>0\) transversely in physical time.

The upper hit occurs at \(H\in[M-\rho,M+2\rho]\). After capture the orbit remains in the tube until \(M-\rho\). At \(t\geq M+\rho\) the speed toward \(u=0\) is at least \(\rho/\epsilon\), and \(2\epsilon/\rho<\rho\). Failure of height to increase to the specified value in finite physical time inside the bounded strip would force the invariant limit set into \(x=0\), hence to the origin. The tube, the negative launch bounds outside the tube and the terminal positive heights exclude this. The upper hit is transverse since \(H>0\).

Write \(H_s^*\) for the hit from \(-1\), \(h_s(r)\) for the one from \(u=1/8,t=r\), and \(H_s^b\) for the one from \(-1/2\). Define \(g_s(r)=H_s^*-h_s(r)\) and \(g_s^b=H_s^*-H_s^b\). For the latter comparison use \(r=-1/2\) below. At the later launch ordinate \(r\), the reference is already in its tube and its \(z\) coordinate exceeds 1.22. The difference \(D_u(r)\) of the two \(u\) positions lies in \((1/2,2)\). The reference stays above the other orbit by scalar order. Up to the first inward crossing of \(u_c=(1.10)^2/2\) by the lower orbit, their positive difference satisfies

\[
 \epsilon(\log D_u)'=[p]_{u_*,u}.
\]

During \([r,r+\rho]\) the secant has magnitude bounded by \(K\). From \(r+\rho\) to \(M-\rho\) it differs from \(p'(U(t))\) by at most \(K\nu\). On \([0.002,0.003]\) the secant is below \(-0.01\), since \(P(1.22)>0.01\) and the tube lies above 1.22. To justify the last statement, the root's distance from \(u=(1.22)^2/2\) exceeds \(\nu\) by the bound on the slope over that interval. On the remainder up to the first crossing of \(u_c\), the elapsed height is at most \(3\rho\) and the secant has magnitude at most \(K\). The difference \(D_c\) at that crossing is less than \(2\exp(-5\cdot10^{-6}/\epsilon)\), using \(4K\rho<5\cdot10^{-6}\). In particular \(D_c<((1.11)^2-(1.10)^2)/2\) by \(e^{-v}<1/v\) and the scale inequality. The reference also lies below \(z=1.11\) at this time. Both traverse \(u_c\) with inward speed between \(c/\epsilon\) and \(K/\epsilon\). Their crossing heights differ by a value between \(\epsilon D_c/K\) and \(\epsilon D_c/c\).

After the respective crossings use \(u\) as the independent variable. Both heights solve \(t_u=-\epsilon/(t-p(u))\) with denominator greater than \(c\). For their positive difference \(d_t\) one has

\[
 \frac{d}{du}\log d_t=
 \frac{\epsilon}{(t_1-p(u))(t_2-p(u))}.
\]

Thus the gap at \(u=0\) lies between \((\epsilon/K)D_c\exp(-\epsilon u_c/c^2)\) and \((\epsilon/c)D_c\). The loss in \(\epsilon\log\) units from this conversion, together with the initial \(\epsilon\log D_u(r)\), has absolute value less than \(6000\epsilon<K\rho\). Here \(\epsilon>10^{-800}\) and \(u_c<1\). The omitted action integrals at the beginning and end of tracking have height lengths \(\rho\) each and root slopes bounded by \(K\). The tracking error is at most \(2K\nu\) and the exit layer costs at most \(3K\rho\). Consequently the corresponding \(g\) satisfies

\[
 \left|-\epsilon\log g-J_P(r)\right|<10K\rho,\qquad
 J_P(r)=\int_{z_m}^{v_P(r)}\frac{P'(z)^2}{z}\,dz ,
\]

where \(z_m=\sqrt{2u_m}\) and \(v_P(r)=\sqrt{2U(r)}\). The identity follows from the change \(t=P(z)\) in \(-\int_r^M p'(U(t))dt\). `falling_action` encloses both endpoints using one-sided sign bisections. It integrates the polynomial terms of \(P'^2/z\) and bounds the remaining \(P'(0)^2\log(v_P/z_m)\) using \(0<\log(v_P/z_m)<1\).

To compare the cap heights, \(P(d)\leq M\leq P(d)+m_\mathrm{lift}\) where \(m_\mathrm{lift}=(10^5\eta)^2/(2\cdot0.1)\). Use \(E'(d)=0\), the bound on \(|p'(d^2/2)|\) from the perturbation alone, and the quadratic upper bound from \(p''\leq-0.1\). The difference of the actual reference hits differs from \(P_+(d)-P_-(d)\) by at most \(m_\mathrm{lift}+3\rho\). The odd value at \(d\) uses the exact identity \(B(d)=\beta\theta\).

For the fixed diagnostic, the rational midpoint of \(P_+(d)\) and \(P_-(d)\) is below \(H_+^*\) and above \(M_-+2\rho\). The plus side returns by the upper trap. In forward direction the minus side reaches \(u=2\) before dropping to \(t=M_-+\rho\), since \(du/d(-t)\geq\rho/\epsilon\) while \(u\leq2\) and \(t\geq M_-+\rho\). For \(u\geq2,t\geq0\) the tail gives \(du/d(-t)\geq10^{-3}(2u)^6/\epsilon\). Its reciprocal has total integral over \(2\leq u<\infty\) less than \(\epsilon\). The orbit escapes before \(t=0\).

## 2. Resolving entrance heights in physical units

The absolute action tolerance in §1 suffices for the two faces in \(d\). A physical height shift has size \(1/S\) in \(t\). To resolve it, factor out the common terminal map on each side.

Use \(b=0.003\) and \(|r|\leq25/S\). A section launch reaches \(z=1.22\) in less than \(1500\epsilon\) in height, since \(p(u)-t>0.0009\) up to that point. The intervening distance in \(u\) is below 1. From then to \(b\) the wall at \(z=1.22\) points right and the wall at 2 points left. The reference is already inside them. The root has \(|U'|\leq100\) and the error equation relative to the root has damping at least \(0.01/\epsilon\). Variation of constants gives \(|u-U(t)|\leq10002\epsilon\) for both arcs after \(r+400000\epsilon\) and up to \(b\). The transient from an initial absolute error at most 2 is bounded by \(2e^{-3000}<2\epsilon\). The recovery term is bounded by \(10^4\epsilon\). Hence

\[
 \left|\epsilon\log D_u(b;r)-\int_r^b p'(U(t))\,dt\right|<10^{13}\epsilon .
\]

The initial gap contributes at most \(2\epsilon\), the launch layer at most \(2K(400000\epsilon)\), and the remainder at most \(10002K\epsilon\) because \(b-r<1\). The section orbits and reference at \(b\) lie in an interval of width at most \(2e^{-5\cdot10^{-6}/\epsilon}\) by the falling-branch contraction on \([0.002,0.003]\).

Let \(K_s(\xi)\) be the upper hit of a solution starting at \(u=\xi,t=b\) in that interval. Order and the tube bound extend the passage estimates to every intermediate launch. Differentiating to a positive spatial cutoff and passing to the transverse exit gives

\[
 \epsilon\log K_s'(\xi)=\epsilon\log(\epsilon/K_s(\xi))
             +\int_b^{K_s(\xi)}p'(u(t;\xi))\,dt .
\]

The factor \(\epsilon/K_s(\xi)\) follows from \(u_t(H)=-H/\epsilon\). To justify the limit from positive cutoff, the last segment uses \(u\) as independent variable and \(t-p(u)>c\). In that chart \(\int_0^{u_c}|p'(u)|du=\int_0^{1.10}|P'(z)|dz<2K\).

Until the first exit at \(u_c\), any pair from the interval stays at distance at most \(D_*=2e^{-2\cdot10^{-6}/\epsilon}\). On the tube their distance can only decrease, and the height range outside the tube before exit is at most \(3\rho\). The two crossing heights differ by at most \(\epsilon D_*/c\). The difference of the integrals in the display up to the earlier crossing is at most \(K D_*\), since the height interval has length below 1. The piece up to the later crossing costs at most \(K\epsilon D_*/c\). In the spatial chart on \([0,u_c]\) the height difference contracts, so the difference of the two integrals there is bounded by \(2K\epsilon^2D_*/c^3\). The prefactor contributes at most \(100\epsilon^2D_*/c\) since both upper hits exceed \(0.01\). The resulting oscillation bound is \(5KD_*\). The inequalities in `fine_section_checks` give \(D_*<\epsilon^2\) by \(e^{-v}\leq6/v^3\).

By the mean value theorem, \(g_s(r)=K_s'(\xi_r)D_u(b;r)\) for some point in the interval. Subtracting the logarithms for two launches gives

\[
 \left|\epsilon\log\frac{g_s(r_1)}{g_s(r_0)}
          +\int_{r_0}^{r_1}p'(U(t))\,dt\right|<3\cdot10^{13}\epsilon .
\]

For each fixed \(d\), let \(v_0>1.22\) be the falling root of \(E(v_0)=0\) and put \(c_d=-E'(v_0)/v_0>0.01\). The perturbation changes the falling root in \(u\) by at most \(100(\|P-E\|+25/S)\), by the lower slope bound. On \(z\geq1.22\), \(|p'(u)-(E'(z)/z)|\leq\|P'-E'\|\). The rational bounds in `fine_section_checks` yield, for \(y_0,y_1\in[-25,1]\),

\[
 \left|S\epsilon\log\frac{g_s(y_1/S)}{g_s(y_0/S)}-c_d(y_1-y_0)\right|<10^{-55}.
\]

At \(H_+^*=H_-^*\) and \(g_+(\Delta/S)=g_-(0)\), use this inequality at \(y_0=\Delta=3/500\) on the plus side and \(y_0=0\) on the minus side. Testing \(y_-=y-\Delta\pm10^{-12}\) gives opposite signs of the gap difference at \(y=4/5,0,-20\). The section hit decreases strictly in the launch height because earlier launches lie to the right of the section when a later launch begins, by the rightward boundary at \(u=1/8\) in this chart. Thus the matching minus ordinate exists uniquely between the two tests. The paired upper hits occur in increasing order as \(y\) descends.

## 3. Compact integrations and the common parameter box

`compact.py` uses the Taylor paths and exact polynomial defects in the included reference copy of `verify_seed.py`. For a reflected sextic \(Q_s\) the completed-square checks imply \(Q_s'\geq-1\) on \(x\geq-1/100\) and \(Q_s'\geq Q_s'(0)>-1/250\) on \(x\geq0\). The symmetric part of the Jacobian of \((y-Q_s(x),-x)\) is \(\operatorname{diag}(-Q_s'(x),0)\). For an approximate path in normalized time over a step \(h\), with exact residual bound \(R\), the Euclidean error is bounded by

\[
 (e_\mathrm{in}+R+h\delta)/(1-\gamma h).
\]

Here \(\delta\) bounds the difference of the fixed field from its sextic on the verified compact domain and \(\gamma=1\), or \(\gamma=1/250\) after checking the broader tube stays in \(x>0\). This follows by integrating \(e'\leq h\gamma e+R+h\delta\) and applying \(e^{h\gamma}\leq(1-h\gamma)^{-1}\). A first-exit argument works because the computed tube stays strictly inside the domain. The normalized-time paths have dyadic coefficients on the \(2^{-256}\) grid. Errors are rounded up. The last step proves strict negative \(x\)-velocity and a sign change in \(x\). Dyadic bisection of the proposed path, the tube and the speed bound enclose the true event time and ordinate. The program checks departure and positivity through every earlier step.

The computation writes \(I_s(h)\) for the quantity denoted by \(\Lambda_s(h)\) in the supplement. Its definition is \(I_s(h)=\int G_s'(x(\tau))\,d\tau\) over the half-arc. The integrator exactly composes the derivative of the sextic with the proposed path and integrates the resulting polynomial. Over a step the derivative error is bounded by \(\|Q_s''\|e_\mathrm{out}+\delta\), using the derivative perturbation as well. The final step includes the product of the event-time error and a bound on \(|G_s'|\). The initial height interval gives one tube for every launch in the interval. Liouville's determinant formula gives

\[
 w_s'(h)=\frac{h}{w_s(h)}e^{-I_s(h)},\qquad
 (w_-^{-1}\circ w_+)'(h)=e^{-(I_+(h)-I_-(h))}
 \quad\text{at }w_+(h)=w_-(h).
\]

The endpoint signs and the nonzero enclosure of \(I_+-I_-\) on each of the five intervals prove a simple fixed point in each. A second fixed point in an interval would force a crossing of opposite orientation, since every zero there has the same nonzero derivative sign for \(D\).

For the finite family, `finite_inner_geometry` proves \(S T>9\), where \(T=(E(1)+E(d))/2\). The wall \(z=1.18\) confines all forward upper launches up to \(S T\). One has \(G_s(x)>8x^6\) on \(2\leq x\leq L/2\). The region \([2,10]\) uses the sextic bound and the physical coefficient error. In \([10,L/2]\), \(E(z)\geq z^6(9/16)(1-1/(4d^2))/6\) from the formula for \(E'\). The code subtracts the possibly negative lower powers and the remaining perturbation. On \(L/2\leq x\leq1.18L\), \(G_s(x)>x^6/100\) by \(E(z)\geq E(1/2)\) and the bound on \(P-E\).

At \(8\leq y\leq9\) the equation for \(v=x^2/2\) is \(dv/d(-y)=y-G_s(x)\). While \(v\geq50\) this is at most \(-v^3/25\), and while \(2\leq v\leq50\) it is at most \(-4v^3\). Their total height decrease is less than \(25/(2\cdot50^2)+1/32<1\). Thus the reference from \(S T\) hits \(y=8\) with \(0<x\leq2\). From any \(q\in[0,2]\) at \(y=8\), the physical wall \(x=3\) and the argument at the origin give a lower-axis hit. Order in the initial \(q\) at \(y=8\) reduces the matching faces to the validated compact arcs from \(q=0,2\). Also from \((L/2,y)\) with \(y\leq1\), the inward spatial equation places the ordinate at \(x=2\) in \([y-1/448,y]\), since \(0<dy/dx=x/(G_s(x)-y)\leq1/(7x^5)\). Compact hit magnitudes from \(x=2\) decrease in the initial ordinate in this range by transversality and order.

For a fixed starting point, determinant variation in the linear coefficient \(a\) gives

\[
 w_s\,\partial_a w_s=-s\int x(\tau)^2
      \exp\!\left(-\int_\tau^{\tau_s} G_s'(x(v))\,dv\right)d\tau .
\]

The same formula holds for starts at fixed \(q>0\). It follows by solving the linearized determinant equation with forcing \(-s x^2\). Thus the half-difference at the inner reference decreases strictly in \(a\), and the shooting faces have the checked strict signs over the full range of the other parameters. Parameter monotonicity also gives the five uniform compact signs. On the three section samples the \(1/448\) ordinate tolerance and the \(10^{-12}\) pairing tolerance allow comparisons at two appropriate compact corners. This gives signs \(-,+,-\) at the paired upper heights of the same matched field.

The three continuous residuals have strictly outward signs on the full rectangular box when ordered as

\[
 \left(\epsilon\log\frac{g_+(\Delta/S)}{g_-(0)},\ H_+^*-H_-^*,\ -D(S T)\right)
\]

in the coordinates \((d,\theta,a)\). Continuity follows from transverse hits and positive gaps. For any positive \(\gamma\), projecting \(p-\gamma R(p)\) onto the box gives a continuous self-map. Brouwer gives a fixed point, the strict face signs put it in the interior, and then \(R(p)=0\). This proves a common real parameter at \(\lambda=10^{-40}\). For fixed \(d,\theta\), the zero \(a_*(d,\theta)\) of \(D(S T)\) is unique in the given \(a\) interval.

## 4. Direct final sign and necessary precision

For either side deform \(E\) to \(P=E+R_s\) by \(P_\gamma=E+\gamma R_s\), \(0\leq\gamma\leq1\). Here \(R_s=P-E\) is the profile perturbation. The checked curvature and falling slope bounds hold throughout the homotopy by convexity. At fixed \(r=-1/2\), differentiate \(J_{P_\gamma}(r)\). Its derivative is

\[
 2\int_{z_m(\gamma)}^{v(\gamma)} \frac{P_\gamma'(z)}{z}R_s'(z)\,dz
       -\frac{P_\gamma'(v(\gamma))}{v(\gamma)}R_s(v(\gamma)).
\]

The peak term vanishes because \(P_\gamma'(z_m)=0\). The falling root has \(|v'(\gamma)|\leq100\mu\) for \(\mu=10^7\eta\), which bounds \(|R_s|\) and \(|R_s'|\). A second differentiation is bounded by \(10^5K\mu^2\). The differentiated integral has interior size at most \(2\mu^2\), zero peak contribution, and falling boundary contribution at most \(200K\mu^2\). In the final term the spatial derivative of \(P_\gamma'/z\) is bounded by \(2K\). Together with \(|v'|\leq100\mu\) this bounds the remaining terms by \((1+300K)\mu^2\). The absolute value of the combined first-order Taylor remainder is at most \(10^5K\mu^2\). At \(\gamma=0\), the difference of the two first derivatives is \(2\eta\int_d^{v_0(r)}(A C'-A'C)dz\) by integration by parts, where \(A=E'/z\) and \(C=(P_+-P_-)/(2\eta)\). `terminal_checks` evaluates this expression and includes both the Taylor remainder and the tolerance from §1. It proves \(g_+^b>g_-^b\). At the common cap match the plus direct hit is therefore lower than the minus direct hit. At the latter upper height \(w_+>S/2=w_-\), so the return sign is positive. The absolute gap bounds put both direct hits above all three section hits and below \(H^*\).

At a common match the five compact heights followed by the three section heights, the minus direct height and the common reference give signs \((+,-,+,-,+,-,+,-,+,0)\). The compact heights lie below the section heights by \(H\in[M-\rho,M+2\rho]\) and the bound on \(S\). Holding all those heights fixed, an upward perturbation of \(a\) gives a negative sign at the last one and retains the strict signs by openness. The perturbation may stay inside the box. At zeros of \(D(h,a)\), \(\partial_aD<0\), so the zero set is a smooth curve. Sard's theorem for its projection to the \(a\) axis gives a nearby regular value where the nine crossings are all simple. Their persistence gives an open set of exact-degree-twelve coefficients containing rational points. This existence argument does not choose a particular rational point.

The actual gap bounds constrain any fixed field in the box using these sections. A common section hit forces \(|H_+^*-H_-^*|\leq\max(g_+,g_-)<\exp(-4\cdot10^{777})\). At a cap match the direct sample and reference differ by less than \(\exp(-3\cdot10^{780})\) in scaled height. Disjoint absolute enclosures of such heights on a uniform binary mesh require more than \(4\cdot10^{777}\) or \(3\cdot10^{780}\) fractional bits respectively, since \(e>2\). This is a bound for that representation of the comparisons.

There is also a necessary condition at the inner reference. `shoulder_obstruction` checks that every listed shoulder height is in \((0.0108,0.0110)\) in scaled units and below \(T\). From any such launch and from \(T\), both arcs on one side enter the strip \(1/2\leq z\leq0.9\) before \(t=0.006\). For \(0.006\leq t\leq0.007\) the speed \(du/d(-t)\) points toward the strip with magnitude at least \(10^{-4}/\epsilon\) outside. The strip is invariant down to \(t=0.002\), with \(p'>0.005\) there. The trajectories contract by \(\exp(-2\cdot10^{-5}/\epsilon)\) on this common interval. Afterward \(x\geq0.1\) holds down to \(y=8\) by the barrier \(G_s(0.1)<8\), and \(z\leq0.9\) also persists. The even leading slope is nonnegative. The verified higher terms together with the linear term evaluated on \(z\geq0.1/L\) give \(p'>-10^{-8}\). The separation at \(y=8\) in physical \(x^2/2\) units is therefore less than \(2L^2\exp(-1.9\cdot10^{-5}/\epsilon)\).

The compact terminal map from \(y=8\) has derivative with respect to \(v=x^2/2\) equal to \(e^{-I_s}/w_s\) by the determinant formula. For \(0\leq x\leq2\) on this section the return magnitudes lie in \((0.03,0.1)\). The endpoint trajectory from \(x=0,y=8\) returns in under 41 time units by the accepted integration. This time bound holds over the full parameter box by validating the arc with an enlarged error of \(10^{-29}\) about the appropriate endpoint sextic. Each larger trajectory takes at most that long to cross the first endpoint's terminal height, since its \(x(y)\) is larger over the common interval. After that height \(\dot x\leq-0.03+0.015\), since \(G_s'(x)>-0.005\) on \(0\leq x\leq3\) and \(G_s(0)=0\). The remaining time is less than 200. Thus the derivative of the terminal map is bounded above by \(e^{1.25}/0.03<200\). The two forward half-return gaps relative to \(T\) are each less than \(\exp(-4\cdot10^{775})\) by the scale checks. Opposite signs at two shoulder heights require \(|D(S T)|<\exp(-4\cdot10^{775})\).

Restrict the determinant formula for \(\partial_a w_s(S T)\) to the compact segment immediately after hitting \(y=8\), where the reference has \(x\geq0.1\). For \(5\cdot10^{-7}\) time units \(x\geq0.05\) since \(|\dot x|<10^5\) on this segment. From there to the hit the time is less than 250 and \(|G_s'|<10^5\). The check \(w_s<0.1\) gives \(-\partial_aD(S T)>\exp(-3\cdot10^7)\). The mean value theorem forces \(|a-a_*(d,\theta)|<\exp(-3\cdot10^{775})\) for the mixed signs at fixed \(d,\theta\). This estimate uses the actual half-returns. The individual signs at the eight shoulder samples still require their finite tracking comparisons.

## 5. A fixed rational field with a sixth hyperbolic cycle

Use \(d=6/5,\ \theta=-10^{-20},\ \lambda=10^{-40}\) and the same exact midpoint \(a=(a_-+a_+)/2\). `fixed_outer.py` checks its exact coefficients against `inputs/outer-variant-field.json`. It also checks the rational upper-height enclosure

\[
 \mathcal I=\big(S(P_+(d)-\rho),\ S(P_+(d)+m_\mathrm{lift}+2\rho)\big),\qquad m_\mathrm{lift}=5\cdot 10^{10}\eta^2,\quad\rho=\lambda^4.
\]

Let \(H_1=S h_+(9/S)\) and \(H_2=S h_+(-1/S)\). Here the argument of \(h_+\) is the scaled height on \(z=1/2\). The earlier launch lies farther right when the later one starts, giving \(H_1<H_2\) by order. Both hits are strictly inside \(\mathcal I\). At \(M_+-\rho\) each backward arc is still in the falling tube at positive \(u\), and from \(M_++\rho\) the remaining height gain is at most \(2\epsilon/\rho<\rho\). Apply \(P_+(d)\le M_+\le P_+(d)+m_\mathrm{lift}\) from §1.

At the rational upper endpoint write \(t_{\rm up}=\sup\mathcal I/S\). Since \(t_{\rm up}\ge M_++2\rho\), the forward equation \(\epsilon\,du/d(-t)=t-p_+(u)\) gives speed at least \(\rho/\epsilon\) while \(u\le2\) and \(t\ge M_++\rho\). Hence the orbit reaches \(u=2\) after a height loss of at most \(2\epsilon/\rho<\rho\). Thereafter \(-p_+(u)>10^{-3}(2u)^6\), so while \(t\ge0\) the remaining loss is bounded by
\[
 \epsilon\int_2^\infty\frac{du}{10^{-3}(2u)^6}=\frac{25}{256}\epsilon<\rho .
\]
Thus \(u\) tends to infinity at positive height. Both \(w_+\) and \(D\) are undefined at this endpoint. The rational check in `../check_fixed_inputs.py` evaluates the JSON endpoints, the cap and tail polynomial bounds, and both loss inequalities.

The code proves \(P_-(d)>P_+(d)+m_\mathrm{lift}+2\rho\). Thus \(x=Ld\) is a forward wall on the minus side for all upper heights in \(\mathcal I\). The lower bound \(G_s(x)>x^6/100\) extends up to \(Ld\) because \(E\) is increasing from \(1/2\) to \(d\) and \(P_s>10^{-3}\) there. The capture argument of §3 gives \(x\leq2\) at \(y=8\) on the minus side. On the plus side the inward portion from \((L/2,9)\) has the same bound at \(y=8\). Its lower hit minus that of the minus side is at most \(B_+(2,8)-B_-(0,8)<-2\cdot10^{-32}\). From the other section launch \((L/2,-1)\) the plus side returns with magnitude more than 1. It reaches \(x=2\) along the inward leg and returns inside the compact wall. The minus magnitude lies in \((0.03,0.1)\) throughout \([H_1,H_2]\) by order and the validated compact arcs. Both half-returns are defined throughout this interval by the returning upper traps. Thus \(D(H_1)<0<D(H_2)\) gives at least one periodic orbit. Every zero in \((H_1,H_2)\) corresponds to an interior plus launch ordinate \(y\in(-1,9)\) by the monotonicity and continuity of the upper section map.

To prove simplicity at such a zero, let \(I_s\) be the integral of the negative of the half-field divergence as in §3. On any half-arc in the scaled variables it has integrand \(p_s'(u)\,|dt|/\epsilon\), since \(|dt/d\tau|=x/S\). On the plus arc from the launch section to the upper hit, use the backward chart of §1. The initial and exit height layers up to the crossing of \(u_c\) contribute at most \(4K\rho\) in \(\epsilon I_+\). While tracking, \(p_+'<0\), including \(p_+'<-0.01\) on \([0.002,0.003]\), so that part contributes at most \(-10^{-5}\). From \(u_c\) to the upper axis the spatial chart contributes at most \(2K\epsilon/c\). On the forward tail from \(L/2\) to \(x=2\), the bound \(|p_+'|\leq K\) and the height change at a zero, \(y+w_+\leq9.1\), give an upper bound \(10K/S\).

The remaining plus arc inside \(0\leq x\leq3\) takes less than 1000 physical time units. While \(x\geq0.01\), \(y\) decreases at rate at least 0.01, and the total available height change is at most 9.1. After entering \(x<0.01\), the orbit cannot cross back through \(x=0.01\) because \(y\) decreases and \(\dot x=y-G_s(0.01)\) on that boundary. In the narrow strip \(|G_s(x)/x|<0.2\) by the checked seed and derivative perturbation bounds, so the polar angular speed is less than \(-0.9\). The angle lies in an interval of width \(\pi\) and the remaining time is less than 5. The same compact time bound applies to the minus arc after reaching \(y=8\), since it then starts inside \(x\leq2\), with \(-w_-\geq-0.1\), and starts at \(x\geq0.01\) by the lower barrier. The compact walls and \(|G_s'|<10^5\) give an absolute contribution to \(\epsilon I_s\) less than \(10^8\epsilon\) on either side.

On the minus side the initial passage from \(x=0\) to \(x=0.01\) is completed while \(t>0.01\). Indeed the launch is above \(t=0.012\), \(t-p(u)\geq0.009\) for \(x\leq0.01\) and \(t\geq0.01\), and the checked scale bound makes the height decrease less than 0.002. Changing to \(u\) there bounds the absolute contribution by \(2K\epsilon/0.009\). From \(x=0.01\) down to \(y=8\), the boundary \(x=0.01\) points right by \(G_s(0.01)<8\), and the cap wall bounds \(z\leq d\). The even slope \(E'/z\) is nonnegative throughout that range. The higher correction terms and the linear term for \(z\geq0.01/L\) give \(p_-'>-10^{-8}\), for a total contribution greater than \(-0.014\cdot10^{-8}\). Combining these bounds gives \(\epsilon I_+<-9\cdot10^{-6}\) and \(\epsilon I_->-10^{-6}\) at every zero in \((H_1,H_2)\). The multiplier is greater than 1 by §3, so the outer crossing is simple and repelling. The same crossing orientation at every zero implies uniqueness between the section hits. The five earlier intervals are verified for the variant too, and they lie strictly below its rational outer enclosure.
