import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.widgets import Slider, Button, CheckButtons
from scipy.integrate import solve_ivp


R_init = 250
L_init = 500
C_init = 20
A_init = 250
omega_init = 250

def rlc_ode(t, y, R, L, C, A, omega):
    dI_dt = (A * np.sin(omega * t) - R * y[0] - y[1]) / L
    dUc_dt = y[0] / C
    return [dI_dt, dUc_dt]

def update(val):

    R = s_R.val
    L = s_L.val / 1000
    C = s_C.val * 1e-6
    A = s_A.val
    omega = s_omega.val

    # Rozwiązanie równań różniczkowych
    t_span = (0, 3 * 2 * np.pi / omega)
    t_eval = np.linspace(t_span[0], t_span[1], 1000)
    y0 = [0, 0]
    sol = solve_ivp(rlc_ode, t_span, y0, t_eval=t_eval, args=(R, L, C, A, omega))

    t = sol.t
    I = sol.y[0]
    U_C = sol.y[1]
    U_R = R * I
    U_L = L * np.gradient(I, t)
    U_src = A * np.sin(omega * t)


    if checkbox.get_status()[0]:
        line_UR.set_ydata(U_R)
        line_UR.set_xdata(t)
        line_UR.set_visible(True)
    else:
        line_UR.set_visible(False)

    if checkbox.get_status()[1]:
        line_UL.set_ydata(U_L)
        line_UL.set_xdata(t)
        line_UL.set_visible(True)
    else:
        line_UL.set_visible(False)

    if checkbox.get_status()[2]:
        line_UC.set_ydata(U_C)
        line_UC.set_xdata(t)
        line_UC.set_visible(True)
    else:
        line_UC.set_visible(False)

    if checkbox.get_status()[3]:
        line_U.set_ydata(U_src)
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

    ax.set_xlim(0, max(t))
    ax.set_ylim(-A * 1.5, A * 1.5)
    ax_I.set_ylim(min(I) * 1.5, max(I) * 1.5)
    fig.canvas.draw_idle()

def reset(event):
    s_R.reset()
    s_L.reset()
    s_C.reset()
    s_A.reset()
    s_omega.reset()
    update(None)


fig, ax = plt.subplots(figsize=(8, 6))
plt.subplots_adjust(left=0.1, bottom=0.4)


line_UR, = ax.plot([], [], label='Napięcie na rezystorze $U_R(t)$', color='y')
line_UL, = ax.plot([], [], label='Napięcie na cewce $U_L(t)$', color='orange')
line_UC, = ax.plot([], [], label='Napięcie na kondensatorze $U_C(t)$', color='red')
line_U, = ax.plot([], [], label='Napięcie źródła $U(t)$', color='cyan', linestyle='--')

ax_I = ax.twinx()
line_I, = ax_I.plot([], [], label='Prąd $I(t)$', color='blue', linestyle='--')


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
s_C = Slider(ax_C, 'C [uF]', 0.1, 500.0, valinit=C_init, valstep=1)
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


update(None)
plt.show()
