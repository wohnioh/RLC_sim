import matplotlib
matplotlib.use('TkAgg')
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, CheckButtons

# Parametry początkowe
R_init = 250
L_init = 500  # mH
C_init = 20  # µF
A_init = 250
omega_init = 250

# Funkcja aktualizująca wykres
def update(val):
    R = s_R.val
    L = s_L.val / 1000
    C = s_C.val * 1e-6
    A = s_A.val
    omega = s_omega.val

    if omega <= 0:
        return


    T = 2 * np.pi / omega
    num_periods = 3  # Zawsze 3 okresy na wykresie
    t = np.linspace(0, num_periods * T, 1000)

    XL = omega * L
    XC = 1 / (omega * C)
    Z = np.sqrt(R**2 + (XL - XC)**2)
    phi = np.arctan((XL - XC) / R)


    I0 = A / Z
    I = I0 * np.sin(omega * t - phi)


    UR = I * R
    UL = I0 * XL * np.cos(omega * t - phi)
    UC = -I0 * XC * np.cos(omega * t - phi)


    if checkbox.get_status()[0]:
        line_UR.set_ydata(UR)
        line_UR.set_xdata(t)
        line_UR.set_visible(True)
    else:
        line_UR.set_visible(False)

    if checkbox.get_status()[1]:
        line_UL.set_ydata(UL)
        line_UL.set_xdata(t)
        line_UL.set_visible(True)
    else:
        line_UL.set_visible(False)

    if checkbox.get_status()[2]:
        line_UC.set_ydata(UC)
        line_UC.set_xdata(t)
        line_UC.set_visible(True)
    else:
        line_UC.set_visible(False)

    if checkbox.get_status()[3]:
        line_U.set_ydata(A * np.sin(omega * t))
        line_U.set_xdata(t)
        line_U.set_visible(True)
    else:
        line_U.set_visible(False)

    if checkbox.get_status()[4]:
        line_I.set_ydata(I)
        line_I.set_xdata(t)
        line_I.set_visible(True)
    else:
        line_I.set_visible(False)


    t_min, t_max = 0, max(t[-1], 1e-9)
    ax.set_xlim(t_min, t_max)

    y_min, y_max = min(-A, UR.min(), UL.min(), UC.min()), max(A, UR.max(), UL.max(), UC.max())
    if abs(y_max - y_min) < 1e-9:
        y_min, y_max = y_min - 1, y_max + 1
    ax.set_ylim(y_min, y_max)

    I_min, I_max = I.min(), I.max()
    if abs(I_max - I_min) < 1e-9:
        I_min, I_max = I_min - 1, I_max + 1
    ax_I.set_ylim(I_min, I_max)

    fig.canvas.draw_idle()

def reset(event):
    s_R.reset()
    s_L.reset()
    s_C.reset()
    s_A.reset()
    s_omega.reset()
    update(None)


T = 2 * np.pi / omega_init
num_periods = 3
t = np.linspace(0, num_periods * T, 1000)

# Obliczenia początkowe
XL = omega_init * (L_init / 1000)
XC = 1 / (omega_init * (C_init * 1e-6))
Z = np.sqrt(R_init**2 + (XL - XC)**2)
phi = np.arctan((XL - XC) / R_init)

I0 = A_init / Z
I = I0 * np.sin(omega_init * t - phi)

UR = I * R_init
UL = I0 * XL * np.cos(omega_init * t - phi)
UC = -I0 * XC * np.cos(omega_init * t - phi)


fig, ax = plt.subplots(figsize=(8, 6))
plt.subplots_adjust(left=0.1, bottom=0.4)


line_UR, = ax.plot(t, UR, label='Napięcie na rezystorze $U_R(t)$', color='y')
line_UL, = ax.plot(t, UL, label='Napięcie na cewce $U_L(t)$', color='orange')
line_UC, = ax.plot(t, UC, label='Napięcie na kondensatorze $U_C(t)$', color='red')
line_U, = ax.plot(t, A_init * np.sin(omega_init * t), label='Napięcie źródła $U(t)$', color='cyan', linestyle='--')


ax_I = ax.twinx()
line_I, = ax_I.plot(t, I, label='Prąd $I(t)$', color='blue', linestyle='--')


ax.set_title('Napięcia i prąd w obwodzie RLC w funkcji czasu')
ax.set_xlabel('Czas $t$ (s)')
ax.set_ylabel('Napięcie (V)')
ax.legend(loc='upper left')
ax.grid()


ax_I.set_ylabel('Prąd (A)')
ax_I.legend(loc='upper right')


rax = plt.axes([0.8, 0.6, 0.15, 0.15])
labels = ['UR', 'UL', 'UC', 'U', 'I']
visibility = [True, True, True, True, True]
checkbox = CheckButtons(rax, labels, visibility)
checkbox.on_clicked(update)


axcolor = 'lightgoldenrodyellow'
ax_R = plt.axes([0.1, 0.3, 0.65, 0.03], facecolor=axcolor)
ax_L = plt.axes([0.1, 0.25, 0.65, 0.03], facecolor=axcolor)
ax_C = plt.axes([0.1, 0.2, 0.65, 0.03], facecolor=axcolor)
ax_A = plt.axes([0.1, 0.15, 0.65, 0.03], facecolor=axcolor)
ax_omega = plt.axes([0.1, 0.1, 0.65, 0.03], facecolor=axcolor)

s_R = Slider(ax_R, 'R [Ohm]', 0.01, 500.0, valinit=R_init, valstep=1)
s_L = Slider(ax_L, 'L [mH]', 10.0, 1000.0, valinit=L_init, valstep=10)
s_C = Slider(ax_C, 'C [µF]', 0.1, 500.0, valinit=C_init, valstep=1)
s_A = Slider(ax_A, 'A [V]', 1.0, 500.0, valinit=A_init, valstep=1)
s_omega = Slider(ax_omega, '$ω$ [rad/s]', 1.0, 500.0, valinit=omega_init, valstep=1)

resetax = plt.axes([0.8, 0.25, 0.1, 0.04])
button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')

s_R.on_changed(update)
s_L.on_changed(update)
s_C.on_changed(update)
s_A.on_changed(update)
s_omega.on_changed(update)
button.on_clicked(reset)

plt.show()
