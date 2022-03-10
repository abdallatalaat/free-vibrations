"""
This script generates the BMD and deformed shape for a simple beam with SDF.
Abdalla Talaat
"""
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
import numpy as np

l1 = 10 #right span
l2 = 15 # left span

EI = 20000
m = 5
k = 3*EI*(l1+l2)/(l1*l2)**2
omega_n = np.sqrt(k/m)

cohesion = 0
zeta = cohesion/(2*omega_n*m)

omega = omega_n * np.sqrt(1-zeta**2)

u_dynamic = 0
u_0 = 3
u_dash_0 = 0.2

u_max = np.sqrt((u_0)**2+((u_dash_0+u_0*zeta*omega_n)/omega)**2)

t_axis3 = np.array([])
u_axis3 = np.array([])


fig, ax = plt.subplots(3,figsize=(12, 7))
plt.subplots_adjust(left=0.5)

text_k = fig.text(0.15, 0.3, "K = {:.2f}\nω = {:.2f}\nT = {:.2f}\nζ = {:.2f}\nu_max = {:.2f}".format(k, omega, 2*np.pi/omega, zeta, u_max), color='red', weight='black', ha='left', size=18)



def u(t, init_u, init_vel, omega, zeta, omega_n): return (np.exp(-1*zeta*omega_n*t))*(init_u*np.cos(omega*t) + ((init_vel+init_u*zeta*omega_n)/omega)*np.sin(omega*t))
def y(x_lins, u_d, l1, l2):
    res = np.array([])
    for x in x_lins:
        if x>= l2: res = np.append(res,(u_d*(l1+l2)/(2*(l1*l2)**2))*((l1*(x**3)/(l1+l2))-(x-l2)**3+x*((l1**3)-l1*(l1+l2)**2)/(l1+l2)))
        else: res = np.append(res,(u_d*(l1+l2)/(2*(l1*l2)**2))*((l1*x**3/(l1+l2))+x*(l1**3-l1*(l1+l2)**2)/(l1+l2)))

    return res

def plot2(color='red', alpha=0.7, lw=2):
    global ax
    ax[2].set_xlabel('time t')
    ax[2].set_ylabel('displacement u(t)')
    ax[2].set_ylim(-6.75, 6.75)

    if len(t_axis3) > 0:
        ax[2].set_xlim(max(t_axis3[-1] - 10, 0), max(10, t_axis3[-1] + 0.5))
        ax[2].plot([max(t_axis3[-1] - 10, 0), max(10, t_axis3[-1] + 0.5)], [0, 0], color='black', lw=0.5)
        ax[2].plot(t_axis3, u_axis3, color=color, lw=lw, alpha=alpha)


def draw_deformed():

    for i in range(len(ax)-1):
        ax[i].set_aspect('equal', adjustable='box')
        ax[i].axis('off')
        ax[i].set_ylim([-6.75, 6.75])

    #deformed beam
    beam_x = np.linspace(0,l1+l2,100)
    beam_y = y(beam_x,u_dynamic,l1,l2)

    #Hinged supports
    support_left_x = np.array([-0.05*(l1+l2), 0, 0.05*(l1+l2)])
    support_right_x = np.array([(l1+l2)*0.95, (l1+l2), (l1+l2)*1.05])
    support_y = np.array([-0.075*(l1+l2), 0, -0.075*(l1+l2)])

    #moment
    moment = k*u_dynamic*l1*l2/(l1+l2)



    #PLOT 0
    ax[0].plot(beam_x, beam_y, color="black", linewidth=3, zorder=1) #beam
    ax[0].scatter(np.array([l2]), np.array([-1 * u_dynamic]), linewidth=6, color="red", zorder=2)
    ax[0].text(l2, -1*u_dynamic+1, "{:.2f}".format(u_dynamic), ha='center', color='red', weight='black')

    #SUPPORTS 0 and 1
    ax[0].plot(support_left_x, support_y, color="black", linewidth=3)
    ax[0].plot(support_right_x, support_y, color="black", linewidth=3)
    ax[1].plot(support_left_x, support_y, color="black", linewidth=3)
    ax[1].plot(support_right_x, support_y, color="black", linewidth=3)


    #PLOT 1
    ax[1].plot(np.array([0, l1+l2]), np.array([0,0]), color="black", linewidth=3, zorder=1)  # beam
    ax[1].plot(np.array([0, l2, l1+l2]), np.array([0, u_dynamic*-1 ,0]), color="red", linewidth=2, zorder=2)  # beam
    ax[1].text(l2, -1 * u_dynamic + 1, "{:.2f}".format(moment), ha='center', color='red', weight='black')

    #PLOT 2
    plot2()


#SLIDERS

ax_EI = plt.axes([0.15, 0.8, 0.3, 0.03])
EI_slider = Slider(
    ax=ax_EI,
    label='Rigidity EI',
    valmin=1000,
    valmax=100000,
    valinit=EI,
)
ax_m = plt.axes([0.15, 0.75, 0.3, 0.03])
m_slider = Slider(
    ax=ax_m,
    label='Mass m',
    valmin=1,
    valmax=50,
    valinit=m,
)
ax_u = plt.axes([0.15, 0.7, 0.3, 0.03])
u_slider = Slider(
    ax=ax_u,
    label='Initial Displacement u_0',
    valmin=0,
    valmax=5,
    valinit=u_0,
)
ax_v = plt.axes([0.15, 0.65, 0.3, 0.03])
v_slider = Slider(
    ax=ax_v,
    label="Initial Velocity u'_0",
    valmin=0,
    valmax=10,
    valinit=u_dash_0,
)
ax_l1 = plt.axes([0.15, 0.6, 0.3, 0.03])
l1_slider = Slider(
    ax=ax_l1,
    label='Right Displacement l1',
    valmin=1,
    valmax=15,
    valinit=l1,
)
ax_l2 = plt.axes([0.15, 0.55, 0.3, 0.03])
l2_slider = Slider(
    ax=ax_l2,
    label="Left Displacement l2",
    valmin=1,
    valmax=25,
    valinit=l2,
)
ax_coh = plt.axes([0.15, 0.50, 0.3, 0.03])
coh_slider = Slider(
    ax=ax_coh,
    label="Cohesion C",
    valmin=0,
    valmax=0.9999*(2*omega_n*m),
    valinit=cohesion,
)

########################

# Update Button
button_ax = plt.axes([0.15, 0.2, 0.1, 0.05])
upd_button = Button(button_ax, 'Update', hovercolor='0.975')

def update_cycle(val):
    global m, l1, l2, EI, u_0, u_dash_0, k, omega, t_axis3, u_axis3, t, ax, fig, text_k, cohesion, omega_n, zeta
    m = m_slider.val
    l1 = l1_slider.val
    l2 = l2_slider.val
    EI = EI_slider.val
    u_0 = u_slider.val
    u_dash_0 = v_slider.val
    cohesion = coh_slider.val

    k = 3 * EI * (l1 + l2) / (l1 * l2) ** 2
    omega_n = np.sqrt(k / m)
    zeta = cohesion / (2 * omega_n * m)
    if zeta > 0.9999: zeta = 0.9999
    omega = omega_n * np.sqrt(1 - zeta ** 2)
    u_max = np.sqrt((u_0) ** 2 + ((u_dash_0 + u_0 * zeta * omega_n) / omega) ** 2)

    coh_slider.valmax = 0.9999*(2*omega_n*m)
    ax_coh.set_xlim(coh_slider.valmin, coh_slider.valmax)

    ax[2].clear()
    plot2(color='blue', alpha=0.3, lw=1.5)
    t_axis3 = np.array([])
    u_axis3 = np.array([])
    t = 0

    fig.texts = []
    text_k = fig.text(0.15, 0.3, "K = {:.2f}\nω = {:.2f}\nT = {:.2f}\nζ = {:.2f}\nu_max = {:.2f}".format(k, omega, 2*np.pi/omega, zeta, u_max), color='red', weight='black', ha='left', size=18)


upd_button.on_clicked(update_cycle)

def handle_close(evt):
    raise SystemExit()
fig.canvas.mpl_connect('close_event', handle_close)


# m_slider.on_changed(update_cycle)
# l1_slider.on_changed(update_cycle)
# l2_slider.on_changed(update_cycle)
# EI_slider.on_changed(update_cycle)
# u_slider.on_changed(update_cycle)
# v_slider.on_changed(update_cycle)

#MAIN LOOP
t = 0
while True:
    u_dynamic = u(t, u_0, u_dash_0, omega, zeta, omega_n)
    t_axis3 = np.append(t_axis3, [t])
    u_axis3 = np.append(u_axis3, [-1*u_dynamic])
    t += 0.07
    ax[0].clear()
    ax[1].clear()
    draw_deformed()

    plt.pause(0.03)




    # remove unnecessary data
    # if len(t_axis3) > 10:
    #     t_axis3 = np.delete(t_axis3, [_ for _ in range(5)])
    #     u_axis3 = np.delete(u_axis3, [_ for _ in range(5)])


