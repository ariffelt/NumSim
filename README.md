# NumSim
This repository contains the code for the tasks in the lecture "Numerical Simulation" at the University of Stuttgart
![image](https://github.com/ariffelt/NumSim/assets/148349447/04baf247-6279-48d6-80a0-5cb7426396ae)
![grafik](https://github.com/ariffelt/NumSim/assets/78861807/d302ff10-c56f-412d-9233-4c28e5d92e9c)
size()=($N_x+2, N_y+2$)

Partitioning:
Always: UIBegin=UJBegin=VIBegin=VJBegin=0
Start the loops always at UIBegin+1, UJBegin+1, VIBegin+1, VJBegin+1
Like this, our 2dArrays always start at 00, the 00 value is the bottom left ghost value.

Parallelization - Ghost layers:
![grafik](https://github.com/ariffelt/NumSim/assets/78861807/fd14b161-1767-4fe4-ba5c-9695e7dec6b5)
