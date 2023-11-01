# necessary for installation, for making cmake & make run:
# incomplete list - order might not work like this, I did a ton of installing, updating, uninstalling,...

conda install tree

conda install cmake

sudo apt install libglu1-mesa-dev

sudo apt install mesa-libGL-devel (?)

conda install -c conda-forge pugixml

sudo apt install libvtk7-dev

#ich musste außerdem in den Ordner /usr/lib einen FakeLink zu libGL.so erstellen, das war bei mir wo anders drin, das geht mit ln -s ...
