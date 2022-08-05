import tkinter as tk
from tkinter import filedialog as fd

root = tk.Tk()


def browsefunc():
    filename = fd.askopenfilename(filetypes=(('csv files', '*.csv'), ("All files", "*.*")))
    ent1.insert(tk.END, filename)


def getPath():
    global path
    path = ent1.get()
    root.destroy()


msg=tk.Label(root, font=40, text="Please choose your csv file")
msg.grid(row=1, column=1)

ent1 = tk.Entry(root, font=40)
ent1.grid(row=2, column=1)

b1 = tk.Button(root,text="Search...",font=40,command=browsefunc)
b1.grid(row=2,column=2)

b2 = tk.Button(root, text= "Start editing", font=40, command=getPath)
b2.grid(row=4, column=1)

root.mainloop()
