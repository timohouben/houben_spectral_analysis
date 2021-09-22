from calc_tc import calc_tc
import numpy as np
from tkinter import *


def long_enough(Smin, Smax, Tmin, Tmax, Lmin, Lmax):
    tc_list = []
    for i in [Smin, Smax]:
        for j in [Tmin, Tmax]:
            for k in [Lmin, Lmax]:
                tc_list.append(calc_tc(k, i, j))
    tc_min = np.min(tc_list)
    tc_max = np.max(tc_list)
    days_min = tc_min * 10
    days_max = tc_max * 10
    return tc_min, tc_max, days_min, days_max


def button_action():
    Smin = float(Smin_e.get())
    Smax = float(Smax_e.get())
    Tmin = float(Tmin_e.get())
    Tmax = float(Tmax_e.get())
    Lmin = float(Lmin_e.get())
    Lmax = float(Lmax_e.get())
    results = long_enough(Smin, Smax, Tmin, Tmax, Lmin, Lmax)
    string = printer(results)
    msg = Message(master, text=string, width=1000)
    msg.config(foreground="white", bg="black", font=("sans", 18, "bold"))
    msg.grid(row=6, column=2, columnspan=3)


def reset_entries():
    Smin_e.delete(0, END)
    Smax_e.delete(0, END)
    Tmin_e.delete(0, END)
    Tmax_e.delete(0, END)
    Lmin_e.delete(0, END)
    Lmax_e.delete(0, END)


results = long_enough(2e-2, 1e-2, 1e-4, 1e-5, 1000, 1000)


def printer(results):
    return str(
        "The characteristic time scale of your aquifer lies between "
        + "{:.1f}".format(results[0])
        + " and "
        + "{:.1f}".format(results[1])
        + " days.\nThe length of your time series should be at least "
        + "{:.1f}".format(results[0] / 10)
        + " days. You would capture the entire frequency response of your aquifer with a time series of "
        + "{:.1f}".format(results[1] / 10)
        + " days."
    )
    # print("Your time series of measurements should be between " + "{:.1f}".format(results[2]) + " and " + "{:.1f}".format(results[3]) + " days long.")


master = Tk()
master.attributes("-fullscreen", True)
master.configure(background="black")
master.grid_rowconfigure(0, weight=100)
master.grid_rowconfigure(1, weight=100)
master.grid_rowconfigure(2, weight=100)
master.grid_rowconfigure(3, weight=100)
master.grid_rowconfigure(4, weight=100)
master.grid_rowconfigure(5, weight=100)
master.grid_rowconfigure(6, weight=100)
master.grid_columnconfigure(0, weight=100)
master.grid_columnconfigure(1, weight=100)
master.grid_columnconfigure(2, weight=100)
master.grid_columnconfigure(3, weight=100)
master.grid_columnconfigure(4, weight=100)
master.grid_columnconfigure(5, weight=100)
master.grid_columnconfigure(6, weight=100)


# row #0
heading = Message(master, text="Is my time series long enough?", width=2000)
heading.config(foreground="white", bg="black", font=("sans", 70, "bold"))
heading.grid(row=0, column=0, columnspan=7)

# heading.config(bg='lightgreen', font=('times', 24, 'italic'))

# row #1
minimum = Message(master, text="Minimum Estimate", width=200)
minimum.config(foreground="white", bg="black", font=("sans", 18, "bold"))
minimum.grid(row=1, column=3)
maximum = Message(master, text="Maximum Estimate", width=200)
maximum.config(foreground="white", bg="black", font=("sans", 18, "bold"))
maximum.grid(row=1, column=4)


# row #2
St = Message(master, text="Aquifer Storativity [-]: (e.g. 1e-2)", width=220)
St.config(foreground="white", bg="black", font=("sans", 18, "bold"))
St.grid(row=2, column=2, sticky="w")
Smin_e = Entry(master)
Smin_e.config(foreground="black", bg="lightblue", font=("sans", 18, "bold"))
Smax_e = Entry(master)
Smax_e.config(foreground="black", bg="lightblue", font=("sans", 18, "bold"))
Smin_e.grid(row=2, column=3)
Smax_e.grid(row=2, column=4)

# row #3
Tt = Message(master, text="Aquifer Transmissivity [m2/s]:", width=300)
Tt.config(foreground="white", bg="black", font=("sans", 18, "bold"))
Tt.grid(row=3, column=2, sticky="w")
Tmin_e = Entry(master)
Tmin_e.config(foreground="black", bg="lightblue", font=("sans", 18, "bold"))
Tmax_e = Entry(master)
Tmax_e.config(foreground="black", bg="lightblue", font=("sans", 18, "bold"))
Tmin_e.grid(row=3, column=3)
Tmax_e.grid(row=3, column=4)

# row #4
Lt = Message(master, text="Aquifer Length [m]:", width=300)
Lt.config(foreground="white", bg="black", font=("sans", 18, "bold"))
Lt.grid(row=4, column=2, sticky="w")
Lmin_e = Entry(master)
Lmin_e.config(foreground="black", bg="lightblue", font=("sans", 18, "bold"))
Lmax_e = Entry(master)
Lmax_e.config(foreground="black", bg="lightblue", font=("sans", 18, "bold"))
Lmin_e.grid(row=4, column=3)
Lmax_e.grid(row=4, column=4)

# row #5
calculate = Button(
    master,
    text="Calculate tc and minimum length of time series!",
    command=button_action,
    width=48,
)
calculate.config(font=("sans", 18, "bold"))
calculate.grid(row=5, column=3, columnspan=3, sticky="w", padx=50)

# row #6
msg = Message(master, text="", width=1000)
msg.config(foreground="black", bg="black", font=("sans", 18, "bold"))
msg.grid(row=6, column=2, columnspan=3)

# row #7
quit = Button(master, text="Quit", fg="green", command=master.destroy)
quit.config(font=("sans", 18, "bold"))
quit.grid(row=7, column=3)
reset = Button(master, text="Reset", command=reset_entries)
reset.config(font=("sans", 18, "bold"))
reset.grid(row=7, column=4)
master.mainloop()
