from Nim import Nim
from transfinite import Ordinal, w
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from tkinter import *

eval_expr = ""  # experession to evaluate using
display_expr = ""  # Global expression string
arith = "ord"  # either 'ord' or 'nim'


def press(key):
    global display_expr
    global eval_expr

    def eval_convert(key_press):  # might be difference in evaluation syntax
        if key_press == "^":
            return "**"

    eval_expr += str(eval_convert(key))
    display_expr += str(key)
    display.set(display_expr)


def equal():
    global display_expr
    try:
        result = str(eval(display_expr.replace("^", "**")))
        display.set(result.replace("**", "^"))
        display_expr = ""
    except:
        display.set("error")
        display_expr = ""


def clear():
    global display_expr
    display_expr = ""
    display.set("")


if __name__ == "__main__":
    root = Tk()
    root.configure(bg="light green")
    root.title("Nimber Calculator")
    root.geometry("270x200")

    display = StringVar()
    entry = Entry(root, textvariable=display)
    entry.grid(columnspan=4, ipadx=70)

    # Display Frame for rendered latex
    tex_frame = Frame(root, relief="sunken", borderwidth=2)
    tex_frame.grid(row=7, column=0, padx=10, pady=10)
    fig = Figure(figsize=(5, 1.5))
    ax = fig.add_subplot()
    canvas = FigureCanvasTkAgg(fig, master=tex_frame)
    canvas.get_tk_widget().grid(row=9, column=0, padx=10, pady=10)

    expr_tex = eval(eval_expr)._repr_latex_() if eval_expr else ""
    try:
        # Render the text in the center of the axes
        ax.text(
            0.5,
            0.5,
            expr_tex,
            horizontalalignment="center",
            verticalalignment="center",
            fontsize=20,
        )  # Adjust font size
    except Exception as e:
        # Show error in the display if LaTeX is bad
        ax.text(
            0.5,
            0.5,
            f"Error: {e}",
            horizontalalignment="center",
            verticalalignment="center",
            fontsize=12,
            color="red",
        )

    # Turn off the axes (we just want the text)
    ax.axis("off")

    # Redraw the canvas
    canvas.draw()
    # 1. Define common styles (your **kwargs)
    number_props = {"fg": "black", "bg": "red", "height": 1, "width": 7}

    operator_props = {
        "fg": "black",
        "bg": "#FF9500",  # A nice orange
        "height": 1,
        "width": 7,
    }

    # 2. Define button layouts in lists
    #    Format: (text, row, col)
    number_layout = [
        ("1", 2, 0),
        ("2", 2, 1),
        ("3", 2, 2),
        ("4", 3, 0),
        ("5", 3, 1),
        ("6", 3, 2),
        ("7", 4, 0),
        ("8", 4, 1),
        ("9", 4, 2),
        ("0", 5, 0),
    ]

    operator_layout = [("+", 2, 3), ("-", 3, 3), ("*", 4, 3), ("/", 5, 3), ("^", 6, 3)]

    # 3. Create buttons by looping through the layouts

    # --- Create Number Buttons ---
    for text, row, col in number_layout:
        # We use 'val=text' in the lambda to "capture" the
        # current value of 'text' in the loop.
        btn = Button(
            root,
            text=text,
            command=lambda val=text: press(val),
            **number_props,  # <-- Unpack the style dictionary
        )
        btn.grid(row=row, column=col, padx=2, pady=2)

    # --- Create Operator Buttons ---
    for text, row, col in operator_layout:
        btn = Button(
            root,
            text=text,
            command=lambda val=text: press(val),
            **operator_props,  # <-- Use the operator styles
        )
        btn.grid(row=row, column=col, padx=2, pady=2)

    # --- Add your special \u03c9 button ---
    # This is now trivial. Just define its properties and add it.
    omega_props = {"fg": "white", "bg": "blue", "height": 1, "width": 7}
    omega_btn = Button(
        root,
        text="\u03c9",  # The unicode character
        command=lambda: press("w"),  # Or 'omega', etc.
        font=("Arial", 12),  # Can add unique properties
        **omega_props,
    )
    omega_btn.grid(row=6, column=0)

    # Other buttons
    eq = Button(root, text="=", fg="black", bg="red", command=equal, height=1, width=7)
    eq.grid(row=5, column=2)
    clr = Button(
        root, text="Clear", fg="black", bg="red", command=clear, height=1, width=7
    )
    clr.grid(row=5, column=1)

    root.mainloop()
