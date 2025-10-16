# welcome to *globr*.

This is a simple, made-from-scratch n-body simulator written in C++ for an ASTR 513 midterm project. It utilizes the computational efficiency of the Barnes-Hut algorithm [(Barnes & Hut, 1986)](https://doi.org/10.1038/324446a0) to create simulations that can support hundreds to $~10^4$ particles.

> (it's possible we could do more, but that hasn't been tested yet. stay tuned!)

I'll probably be back in the future to play around with this more and add features + updates, but for now, enjoy and feel free to explore!


## how do I use *globr*?

**devnote >>** *globr* has only been built and run on M3 Mac and Ubuntu. It requires at minimum C++11, but uses no libraries outside of the standards you'd expect with a base installation.

1. compile the code! 
    *globr* has a Makefile ready to go, so compilation should be as easy as
    ```
    cd path/to/globr
    mkdir data
    cd src
    make clean; make; ls
    ```
    after you list your directory you should see an executable called `globr`. if so, congrats! compilation worked and you're ready to make cool clusters.
    
    **devnote >>** you might get some warnings about deprecation if you're on a Mac depending on your version of C/C++; ignore those.

2. choose your starting conditions.
    *globr* has a simple set of command line flags that you can set when running the main executable:
    
    - `-N`: number of particles (REQUIRED)
    - `--size`: starting size of physical simulation space, in parsecs (REQUIRED)
    - `--step`: size of timestep in years (default: 1)
    - `--nstep`: number of timesteps (default: 5000)
    - `--freq`: how often data is output, in timesteps (default: 5)
    - `--theta`: barnes-hut criterion (default: 0.5)
    - `--run`: name of your simulation run (and data directory) (REQUIRED)
    - `--init`: initial conditions file (REQUIRED)

    so if we wanted to run 
    - a cluster of $10^4$ stars, 
    - made from the Salpeter initial condtions (salpeter.txt) 
    - for 5000 timesteps,
    - with a step size of 100 years (or 500,000 years) 
    - and output data files every 10 timesteps, 
    then we would run globr from `globr/src` as follows:

    ```
    ./globr -N 10000 --size 100 --step 15 --nstep 5000 --freq 10 --run salpeter --init salpeter.txt &
    ```

3. run *globr* and wait...
    Put together what you've learned in the previous steps and make some clusters! 

    **devnote.** If you try to run *globr* and immediately get a segfault error, try making your simulation size larger. I haven't added dynamic rescaling when the first tree is being constructed yet; if a particle is outside of that domain, the program will crash.

4. time for pretty pictures!
    So now you've made a cluster. What next? Visualization, of course! In `globr/viz` there's a simple Jupyter Notebook called `viz.ipynb`. Before you get started, create a new directory called `frames`:

    ```
    cd globr/viz
    mkdir frames
    ```

    Now you're ready to go. In the notebook, there's two variables:
    `datadir` and `rad_multiplier`.

    `datadir` should be set to whatever you named your run (i.e. `--run`) when creating your simulation; that will find the output data.

    `rad_multiplier` controls how zoomed in your final visualization will be. The default value is 2; make it lower (~0.5) for more dense or smaller clusters.

    Run the cells! You'll need to have `pyvista` installed on your kernel, and also will need `ffmpeg` installed on your laptop or computing tool of choice. Both of these are publicly available and pretty easy to acquire and install on most operating systems through common command line tools (`pip`, `brew`, etc.)

    Now that you have images (so. many.), go to the new directory (`globr/viz/frames/datadir`) and run the following:

    `ffmpeg -framerate 20 -i frame_%05d.png -c:v libx264 -pix_fmt yuv420p NAME.mp4`

    where **NAME** is the filename you want for the visualization. make sure you have `ffmpeg` installed by running `ffmpeg -version` in your terminal beforehand.

    Enjoy your new cool movie of your cluster(s)!

5. do science? 
    You now know how to use and operate *globr*, congrats! There are a variety of initial conditions files in `globr/init` that can be used to 
    - test scaling (`globr/init/scaling_N*.txt`); 
    - explore different initial mass functions, or IMFs (`globr/init/*`; Salpeter, Kroppa, top-heavy, bottom-heavy...); 
    - watch two clusters collide (`globr/init/initialConditions_2cluster*.txt`); 
    - or even just look at a random cluster (`globr/init/initialConditions_10*.txt`).

    Enjoy!



## questions? need help?

Questions, comments, and suggestions are more than welcome. Reach out to the developer (hi, I'm Logan) at lawhite [at] arizona [dot] edu.
