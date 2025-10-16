import matplotlib.pyplot as plt
import matplotlib.animation as anim
import numpy as np

def make_gif(sol_x, filename = "n body sim", frames = 100, xmax = None, ymax = None): 
    """
    Makes a gif out of the input x coordinate history array. This only makes a plot on the 2D x-y plane. 
    sol_x     - the array containing x coordinates, formatted as (Nsteps, Nparticles, 3). Units are in m. Can change later. 
    filename  - the name of the file of which we should save this plot to. No .gif needed!
    frames    - the number of frames you want to include. This always starts from 0, so if you want to do something else, 
                you'll have to truncate sol_x yourself.
    xlim,ylim - limits on x and y axis. Code knows how to determine it now. Code grew up. 
    """
    # set xlim and ylim first
    if xmax is None:
        maxx = np.max(sol_x[:, :, 0])
        mag = np.floor(np.log10(maxx)) # this tells us how many orders of magnitude we need at least
        top = np.ceil(maxx/(10**(mag-1))) # this is the number * orders of magnitude 
        xmax = top * 10**(mag-1)
        
    if ymax is None:
        maxy = np.max(sol_x[:, :, 1])
        mag = np.floor(np.log10(maxy)) # this tells us how many orders of magnitude we need at least
        top = np.ceil(maxy/(10**(mag-1))) # this is the number * orders of magnitude 
        ymax = top * 10**(mag-1)
    
    fig, ax = plt.subplots(1, 1)
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
    
    ax.set(xlim=[-xmax, xmax], ylim=[-ymax, ymax], xlabel='x offset [m]', ylabel='y offset [m]')
    
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

bb_dict = {}
with open("proj1/bbcolors.txt") as f:
    for line in f:
       (key, val) = line.split()
       bb_dict[int(key)] = val

bbtemp = np.asarray(list(bb_dict.keys()))

def make_gif_bb(sol_x, filename = "n body sim", mass = None, temp = None, frames = 100, xmax = None, ymax = None): 
    """
    Makes a gif out of the input x coordinate history array. This only makes a plot on the 2D x-y plane. 
    The dots are now colored by their blackbody temperatures because we are astronomers.
    
    sol_x     - the array containing x coordinates, formatted as (Nsteps, Nparticles, 3). Units are in m. Can change later. 
    filename  - the name of the file of which we should save this plot to. No .gif needed!
    mass      - mass of the stars, used to calculate their temperature if temperature is not given (if temp is given, 
                that one always take the highest priority).
    temp      - temperature of the stars! If no value is given (None), then temperature is calculated by equating 
                    L = sigma*T^4*4*pi*R^2 = M(Msun)**3.5 (units of 1 Lsun) where R propto M**0.7.
                Then T is found through this. Finally T is passed through a splitting thing to find the closest color in the bb_colors.txt.
    frames    - the number of frames you want to include. This always starts from 0, so if you want to do something else, 
                you'll have to truncate sol_x yourself.
    xlim,ylim - limits on x and y axis. This code is stupid and it does not know how to get a good limit. Please help this code.  
    """
    # Check that mass or temp is given so we can color it. 
    if (mass is None) and (temp is None):
        print("You need to give me at least one of mass or temperature to work with. Please ; - ;\"...")
        return None
    elif temp is None:
        # Calculate temperature given mass
        temp = [bb_dict[find_T_bb(m)] for m in mass]

    # set xlim and ylim first
    if xmax is None:
        maxx = np.max(sol_x[:, :, 0])
        mag = np.floor(np.log10(maxx)) # this tells us how many orders of magnitude we need at least
        top = np.ceil(maxx/(10**(mag-1))) # this is the number * orders of magnitude 
        xmax = top * 10**(mag-1)
        
    if ymax is None:
        maxy = np.max(sol_x[:, :, 1])
        mag = np.floor(np.log10(maxy)) # this tells us how many orders of magnitude we need at least
        top = np.ceil(maxy/(10**(mag-1))) # this is the number * orders of magnitude 
        ymax = top * 10**(mag-1)
    
    plt.style.use('dark_background')
    fig, ax = plt.subplots()
    lines = []
    scats = []
    
    for i in range(sol_x.shape[1]):
        line = ax.plot(
            sol_x[0, i, 0],
            sol_x[0, i, 1], 
            color = temp[i], 
            alpha = 0.3
        )[0]
        lines.append(line)
    
    for i in range(sol_x.shape[1]):
        scat = ax.scatter(
            sol_x[0, i, 0],
            sol_x[0, i, 1], 
            color = temp[i]
        )
        scats.append(scat)
    
    ax.set(xlim=[-xmax, xmax], ylim=[-ymax, ymax], xlabel='x offset [m]', ylabel='y offset [m]')
    
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
    plt.style.use('default')

def find_T_bb(M):
    # M in units of solar mass.
    # T**4 = cM**2.1 in solar: find c to be
    c = 5778**4
    T4 = c*M**2.1
    T = T4**(0.25)
    loc = np.asarray(np.where(bbtemp < T)[0])[-1]
    # if we are out of range, loc = len(bbtemp)-1, then we just use the last entry:
    if loc == len(bbtemp) - 1:
        return bbtemp[-1]
    # Determine which one is closer
    if abs(bbtemp[loc] - T) < abs(bbtemp[loc+1] - T): 
        # Closer to loc
        temp = bbtemp[loc]
    else: 
        # closer to the previous one
        temp = bbtemp[loc+1]
        
    return temp