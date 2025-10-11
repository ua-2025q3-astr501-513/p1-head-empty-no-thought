import matplotlib.pyplot as plt
import matplotlib.animation as anim

def make_gif(sol_x, filename = "n body sim", frames = 1000, xmax = 1e9, ymax = 1e9): 
    """
    Makes a gif out of the input x coordinate history array. This only makes a plot on the 2D x-y plane. 
    sol_x     - the array containing x coordinates, formatted as (Nsteps, Nparticles, 3). Units are in km. Can change later. 
    filename  - the name of the file of which we should save this plot to. No .gif needed!
    frames    - the number of frames you want to include. This always starts from 0, so if you want to do something else, 
                you'll have to truncate sol_x yourself.
    xlim,ylim - limits on x and y axis. This code is stupid and it does not know how to get a good limit. Please help this code.  
    """
    
    fig, ax = plt.subplots()
    lines = []
    scats = []
    
    for i in range(sol_x.shape[1]):
        line = ax.plot(
            sol_x[0, i, 0],
            sol_x[0, i, 1]
        )[0]
        lines.append(line)
    
    for i in range(sol_x.shape[1]):
        scat = ax.scatter(
            sol_x[0, i, 0],
            sol_x[0, i, 1]
        )
        scats.append(scat)
    
    ax.set(xlim=[-xmax, xmax], ylim=[-ymax, ymax], xlabel='x offset [km]', ylabel='y offset [km]')
    
    def update(frame):
        # update the scatter plot:
        for i in range(len(scats)):
            scats[i].set_offsets([sol_x[frame, i, 0], sol_x[frame, i, 1]])
        # update the line plot:
        for i in range(len(lines)):
            lines[i].set_xdata(sol_x[:frame, i, 0])
            lines[i].set_ydata(sol_x[:frame, i, 1])
        return (scats, lines)
    
    ani = anim.FuncAnimation(fig=fig, func=update, frames=frames, interval=30)
    ani.save(f"{filename}.gif")