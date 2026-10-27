# Rzut ukośny bez oporu i z oporem powietrza

*Raport obliczeniowy: porównanie rozwiązania analitycznego z symulacją numeryczną.*

**Dane początkowe:** $v_0=37\,\mathrm{m/s}$, $\alpha=75^\circ$.  
**Układ SI:** położenie w metrach, czas w sekundach, prędkość w metrach na sekundę, przyspieszenie w metrach na sekundę do kwadratu.

## Spis treści

1. [Cel i założenia](#1-cel-i-założenia)
2. [Rzut ukośny bez oporu](#2-rzut-ukośny-bez-oporu)
3. [Rzut ukośny z oporem powietrza](#3-rzut-ukośny-z-oporem-powietrza)
4. [Wyniki i porównanie](#4-wyniki-i-porównanie)
5. [Wykresy](#5-wykresy)
6. [Kod i odtworzenie wyników](#6-kod-i-odtworzenie-wyników)
7. [Wnioski](#7-wnioski)

---

## 1. Cel i założenia

Rzut ukośny jest ruchem w dwóch wymiarach. W polu grawitacyjnym składowa pozioma i pionowa zachowują się odmiennie: bez oporu pozioma prędkość jest stała, a pionowa zmienia się pod wpływem grawitacji. Opór powietrza sprzęga te składowe, ponieważ jego kierunek jest przeciwny do całkowitej prędkości.

W obu modelach przyjmujemy:

- ciało startuje z poziomu podłoża, czyli $x_0=0\,\mathrm{m}$ i $y_0=0\,\mathrm{m}$;
- oś $x$ jest pozioma, a oś $y$ skierowana pionowo w górę;
- przyspieszenie grawitacyjne ma wartość $g=9{,}81\,\mathrm{m/s^2}$;
- Ziemia jest płaska na rozpatrywanym obszarze, a $g$ nie zależy od położenia;
- pomijamy wiatr, wypór powietrza, obrót ciała i zmianę gęstości atmosfery.

Do obliczeń z oporem przyjmujemy modelową kulę o masie $m=0{,}145\,\mathrm{kg}$, promieniu $r=0{,}037\,\mathrm{m}$, gęstości powietrza $\rho=1{,}225\,\mathrm{kg/m^3}$ oraz współczynniku oporu $C_d=0{,}47$. Są to jawne założenia przykładu: sam zestaw $v_0$ i $\alpha$ nie wystarcza do wyznaczenia ruchu z oporem. Model nie uwzględnia zmian $C_d$ związanych z liczbą Reynoldsa.

Pole przekroju kuli wynosi

$$
A=\pi r^2=\pi(0{,}037\,\mathrm{m})^2\approx 4{,}30\cdot10^{-3}\,\mathrm{m^2}.
$$

| Symbol | Znaczenie | Jednostka SI |
| :--- | :--- | :---: |
| $x,y$ | współrzędne położenia | $\mathrm{m}$ |
| $v_x,v_y,v$ | składowe prędkości i jej wartość | $\mathrm{m/s}$ |
| $a_x,a_y$ | składowe przyspieszenia | $\mathrm{m/s^2}$ |
| $t$ | czas | $\mathrm{s}$ |
| $m$ | masa ciała | $\mathrm{kg}$ |
| $g$ | przyspieszenie grawitacyjne | $\mathrm{m/s^2}$ |
| $\rho$ | gęstość powietrza | $\mathrm{kg/m^3}$ |
| $C_d$ | bezwymiarowy współczynnik oporu | $1$ |
| $A$ | pole przekroju poprzecznego | $\mathrm{m^2}$ |

## 2. Rzut ukośny bez oporu

### 2.1. Rozkład prędkości początkowej

Prędkość początkową rozkładamy na składowe wzdłuż osi:

$$
v_{0x}=v_0\cos\alpha,
\qquad
v_{0y}=v_0\sin\alpha.
$$

Dla podanych danych:

$$
v_{0x}=37\cos75^\circ\approx 9{,}576\,\mathrm{m/s},
\qquad
v_{0y}=37\sin75^\circ\approx 35{,}739\,\mathrm{m/s}.
$$

### 2.2. Równania ruchu i ich wyprowadzenie

Bez oporu jedyną siłą jest ciężar $\vec F_g=(0,-mg)$. Z drugiej zasady Newtona $m\vec a=\vec F$ wynika:

$$
a_x=0,
\qquad
a_y=-g.
$$

Korzystając z $a_x=dv_x/dt$ oraz $a_y=dv_y/dt$ i całkując po czasie z warunkami $v_x(0)=v_{0x}$, $v_y(0)=v_{0y}$, otrzymujemy:

$$
v_x(t)=v_{0x},
\qquad
v_y(t)=v_{0y}-gt.
$$

Ponieważ $v_x=dx/dt$ i $v_y=dy/dt$, ponowne całkowanie z $x(0)=y(0)=0$ daje równania toru parametrycznego:

$$
x(t)=v_0\cos\alpha\,t,
\qquad
y(t)=v_0\sin\alpha\,t-\frac{1}{2}gt^2.
$$

Eliminujemy parametr $t$ z pierwszego równania: $t=x/(v_0\cos\alpha)$. Po podstawieniu do równania pionowego otrzymujemy parabolę:

$$
y(x)=x\tan\alpha-\frac{g x^2}{2v_0^2\cos^2\alpha}.
$$

### 2.3. Czas lotu, wysokość maksymalna i zasięg

Chwila osiągnięcia wierzchołka toru wynika z warunku $v_y(t_H)=0$:

$$
t_H=\frac{v_0\sin\alpha}{g}\approx 3{,}643\,\mathrm{s}.
$$

Podstawiając ten czas do $y(t)$, uzyskujemy wysokość maksymalną:

$$
H=\frac{v_0^2\sin^2\alpha}{2g}\approx 65{,}10\,\mathrm{m}.
$$

Poza rozwiązaniem $t=0$, równanie $y(t)=0$ daje czas powrotu do poziomu startu:

$$
t_f=\frac{2v_0\sin\alpha}{g}\approx 7{,}286\,\mathrm{s}.
$$

Zasięg poziomy to $R=x(t_f)$:

$$
R=\frac{v_0^2\sin(2\alpha)}{g}\approx 69{,}78\,\mathrm{m}.
$$

Wartość prędkości w dowolnej chwili wynosi $v(t)=\sqrt{v_x^2+v_y^2}$. Bez oporu składowa pozioma jest stała, zaś pionowa liniowo maleje; w najwyższym punkcie $v_y=0$, lecz $v_x$ nadal jest niezerowe.

## 3. Rzut ukośny z oporem powietrza

### 3.1. Kwadratowe prawo oporu

Dla prędkości typowych dla tego przykładu przyjmujemy opór proporcjonalny do kwadratu szybkości. Wektor siły oporu jest przeciwny do wektora prędkości:

$$
\vec F_d=-\frac{1}{2}\rho C_d A\,|\vec v|\,\vec v,
\qquad
|\vec v|=\sqrt{v_x^2+v_y^2}.
$$

Zapis wektorowy automatycznie nadaje oporowi właściwy zwrot w obu osiach. Z drugiej zasady Newtona, po podzieleniu przez $m$, otrzymujemy:

$$
\frac{dv_x}{dt}=-c|\vec v|v_x,
\qquad
\frac{dv_y}{dt}=-g-c|\vec v|v_y,
\qquad
c=\frac{\rho C_d A}{2m}.
$$

Dla przyjętej kuli współczynnik wynosi

$$
c=\frac{(1{,}225\,\mathrm{kg/m^3})(0{,}47)(4{,}30\cdot10^{-3}\,\mathrm{m^2})}{2(0{,}145\,\mathrm{kg})}
\approx 8{,}54\cdot10^{-3}\,\mathrm{m^{-1}}.
$$

Jednostką $c$ jest $\mathrm{m^{-1}}$, dzięki czemu $c|\vec v|v$ ma jednostkę przyspieszenia $\mathrm{m/s^2}$.

### 3.2. Układ równań pierwszego rzędu

Po dołączeniu definicji prędkości dostajemy układ czterech równań:

$$
\begin{aligned}
\frac{dx}{dt}&=v_x, & \frac{dy}{dt}&=v_y, \\
\frac{dv_x}{dt}&=-c\sqrt{v_x^2+v_y^2}\,v_x,
& \frac{dv_y}{dt}&=-g-c\sqrt{v_x^2+v_y^2}\,v_y.
\end{aligned}
$$

Warunki początkowe to $x(0)=y(0)=0$, $v_x(0)=v_0\cos\alpha$ oraz $v_y(0)=v_0\sin\alpha$. Ze względu na zależność oporu od całkowitej szybkości równania są sprzężone i nieliniowe. Nie mają prostego rozwiązania analitycznego w postaci elementarnych funkcji; dlatego rozwiązujemy je numerycznie metodą Rungego–Kutty czwartego rzędu (RK4).

W każdym kroku RK4 szacuje pochodne czterokrotnie: na początku kroku, w dwóch punktach pośrednich i na końcu. Ważona kombinacja tych pochodnych daje przybliżony stan po czasie $\Delta t$. W kodzie przyjmujemy $\Delta t=0{,}005\,\mathrm{s}$, a lot kończymy w chwili przecięcia poziomu $y=0$ po starcie.

### 3.3. Energia mechaniczna

W modelu bez oporu energia mechaniczna jest zachowana. Opór wykonuje ujemną pracę, zamieniając część energii mechanicznej na energię wewnętrzną powietrza i ciała. Moc siły oporu jest równa:

$$
P_d=\vec F_d\cdot\vec v=-\frac{1}{2}\rho C_d A|\vec v|^3\leq 0.
$$

Zatem energia kinetyczna i potencjalna nie sumują się do stałej podczas lotu z oporem.

## 4. Wyniki i porównanie

Wyniki modelu bez oporu są analityczne. Wartości dla oporu uzyskuje się numerycznie z określonych wyżej parametrów kuli. Skrypt drukuje dokładne wartości symulacji w konsoli, co pozwala je zweryfikować po zmianie parametrów lub kroku czasowego.

| Wielkość | Bez oporu | Z oporem kwadratowym |
| :--- | ---: | ---: |
| Czas lotu | $7{,}286\,\mathrm{s}$ | wynik symulacji |
| Wysokość maksymalna | $65{,}10\,\mathrm{m}$ | wynik symulacji |
| Zasięg | $69{,}78\,\mathrm{m}$ | wynik symulacji |
| Prędkość pozioma | stała: $9{,}576\,\mathrm{m/s}$ | maleje wskutek oporu |
| Energia mechaniczna | zachowana | maleje |

Opór zmniejsza składowe prędkości i skraca zasięg oraz wysokość maksymalną w porównaniu z lotem w próżni. W przeciwieństwie do rzutu bez oporu, tor z oporem nie jest symetryczny względem wierzchołka: podczas opadania ciało ma mniejszą składową poziomą niż podczas wznoszenia na tej samej wysokości.

## 5. Wykresy

Wykres toru porównuje parabolę w próżni z numerycznie obliczonym torem dla kuli w powietrzu.

![Wykres toru rzutu ukośnego bez oporu i z oporem powietrza](trajectory.png)

*Rysunek 1. Trajektorie $y(x)$; obie osie są w metrach. Punktami oznaczono położenie najwyższe.*

Składowe prędkości i wysokość pokazują, jak opór wpływa na ruch w czasie.

![Prędkość pozioma, pionowa i wysokość w funkcji czasu](time_series.png)

*Rysunek 2. Przebiegi w czasie. Bez oporu $v_x$ jest stałe, a $v_y$ zmienia się liniowo.*

Bilans energii pokazuje różnicę między zachowaniem energii w próżni a jej rozpraszaniem przez opór.

![Energia kinetyczna, potencjalna i mechaniczna w funkcji czasu](energy.png)

*Rysunek 3. W modelu bez oporu $E_k+E_p$ pozostaje stałe. Przy oporze suma maleje.*

## 6. Kod i odtworzenie wyników

Kompletny skrypt znajduje się w pliku [`wykresy.py`](wykresy.py). Wymaga Pythona 3 oraz bibliotek NumPy i Matplotlib. W terminalu VS Code, w tym katalogu, można uruchomić:

```bash
python -m pip install numpy matplotlib
python wykresy.py
```

Skrypt zapisuje trzy pliki PNG w tym samym katalogu co raport oraz wypisuje porównanie liczbowe. Poniżej przedstawiono kod używany do obliczeń i rysowania:

```python
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


V0 = 37.0
ANGLE_DEG = 75.0
G = 9.81
MASS = 0.145
RADIUS = 0.037
RHO = 1.225
CD = 0.47
DT = 0.005
OUTPUT_DIR = Path(__file__).resolve().parent
AREA = np.pi * RADIUS**2
DRAG_COEFFICIENT = RHO * CD * AREA / (2 * MASS)
ANGLE = np.deg2rad(ANGLE_DEG)
VX0 = V0 * np.cos(ANGLE)
VY0 = V0 * np.sin(ANGLE)


def derivatives(state):
	x, y, vx, vy = state
	speed = np.hypot(vx, vy)
	return np.array([
		vx,
		vy,
		-DRAG_COEFFICIENT * speed * vx,
		-G - DRAG_COEFFICIENT * speed * vy,
	])


def rk4_step(state, dt):
	k1 = derivatives(state)
	k2 = derivatives(state + 0.5 * dt * k1)
	k3 = derivatives(state + 0.5 * dt * k2)
	k4 = derivatives(state + dt * k3)
	return state + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4)


def simulate_with_drag():
	times = [0.0]
	states = [np.array([0.0, 0.0, VX0, VY0])]
	while states[-1][1] >= 0.0 or len(states) == 1:
		next_state = rk4_step(states[-1], DT)
		next_time = times[-1] + DT
		if next_state[1] < 0.0:
			fraction = states[-1][1] / (states[-1][1] - next_state[1])
			next_state = states[-1] + fraction * (next_state - states[-1])
			next_time = times[-1] + fraction * DT
		times.append(next_time)
		states.append(next_state)
		if next_state[1] == 0.0:
			break
	return np.array(times), np.vstack(states)


def main():
	flight_time = 2 * VY0 / G
	time_vacuum = np.linspace(0.0, flight_time, 1200)
	x_vacuum = VX0 * time_vacuum
	y_vacuum = VY0 * time_vacuum - 0.5 * G * time_vacuum**2
	vx_vacuum = np.full_like(time_vacuum, VX0)
	vy_vacuum = VY0 - G * time_vacuum
	time_drag, state_drag = simulate_with_drag()
	x_drag, y_drag, vx_drag, vy_drag = state_drag.T

	fig, ax = plt.subplots(figsize=(9, 5.5), constrained_layout=True)
	ax.plot(x_vacuum, y_vacuum, label="Bez oporu", linewidth=2.2)
	ax.plot(x_drag, y_drag, label="Z oporem", linewidth=2.2)
	ax.scatter([x_vacuum.max() / 2], [y_vacuum.max()], s=35, zorder=3)
	ax.scatter([x_drag[y_drag.argmax()]], [y_drag.max()], s=35, zorder=3)
	ax.set(title="Tor rzutu ukośnego", xlabel="x [m]", ylabel="y [m]")
	ax.grid(True, alpha=0.3)
	ax.legend()
	fig.savefig(OUTPUT_DIR / "trajectory.png", dpi=180)
	plt.close(fig)

	fig, axes = plt.subplots(3, 1, figsize=(9, 8), sharex=True,
							 constrained_layout=True)
	for ax, vacuum, drag, label in zip(
		axes,
		(vx_vacuum, vy_vacuum, y_vacuum),
		(vx_drag, vy_drag, y_drag),
		("vx [m/s]", "vy [m/s]", "y [m]"),
	):
		ax.plot(time_vacuum, vacuum, label="Bez oporu", linewidth=2)
		ax.plot(time_drag, drag, label="Z oporem", linewidth=2)
		ax.set_ylabel(label)
		ax.grid(True, alpha=0.3)
	axes[0].legend()
	axes[-1].set_xlabel("t [s]")
	fig.savefig(OUTPUT_DIR / "time_series.png", dpi=180)
	plt.close(fig)

	vacuum_kinetic = 0.5 * MASS * (vx_vacuum**2 + vy_vacuum**2)
	vacuum_potential = MASS * G * y_vacuum
	drag_kinetic = 0.5 * MASS * (vx_drag**2 + vy_drag**2)
	drag_potential = MASS * G * y_drag
	fig, ax = plt.subplots(figsize=(9, 5.5), constrained_layout=True)
	ax.plot(time_vacuum, vacuum_kinetic, label="Ek, bez oporu")
	ax.plot(time_vacuum, vacuum_potential, label="Ep, bez oporu")
	ax.plot(time_vacuum, vacuum_kinetic + vacuum_potential,
			label="E, bez oporu", linewidth=2.4)
	ax.plot(time_drag, drag_kinetic + drag_potential,
			label="E, z oporem", linewidth=2.4)
	ax.set(title="Energia mechaniczna", xlabel="t [s]", ylabel="E [J]")
	ax.grid(True, alpha=0.3)
	ax.legend()
	fig.savefig(OUTPUT_DIR / "energy.png", dpi=180)
	plt.close(fig)

	print(f"Pole przekroju: {AREA:.6f} m^2")
	print(f"Współczynnik oporu c: {DRAG_COEFFICIENT:.6f} 1/m")
	print(f"Bez oporu: czas = {flight_time:.3f} s, "
		  f"wysokość = {y_vacuum.max():.3f} m, zasięg = {x_vacuum[-1]:.3f} m")
	print(f"Z oporem: czas = {time_drag[-1]:.3f} s, "
		  f"wysokość = {y_drag.max():.3f} m, zasięg = {x_drag[-1]:.3f} m")


if __name__ == "__main__":
	main()
```

## 7. Wnioski

1. Bez oporu rzut ukośny jest superpozycją ruchu jednostajnego w poziomie i jednostajnie przyspieszonego w pionie, a jego tor jest parabolą.
2. Dla $v_0=37\,\mathrm{m/s}$ i $\alpha=75^\circ$ model próżni daje wysokość maksymalną około $65{,}10\,\mathrm{m}$, czas lotu około $7{,}286\,\mathrm{s}$ i zasięg około $69{,}78\,\mathrm{m}$.
3. Przy oporze kwadratowym siła zależy od szybkości oraz parametrów ciała i powietrza. Składowe ruchu są sprzężone, więc trajektorię wyznaczono metodą numeryczną.
4. Opór zmniejsza zasięg i wysokość oraz powoduje spadek energii mechanicznej. Wyniki dla oporu należy interpretować w kontekście przyjętych parametrów kuli i modelu atmosfery.
