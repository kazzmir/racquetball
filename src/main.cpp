#include <allegro5/allegro.h>
#include <iostream>

#include "pointer.h"

namespace Physics{

class Vector{
public:
    Vector(double x, double y, double z):
    x(x), y(y), z(z){
    }

    Vector operator+(const Vector & what){
        return Vector(x + what.x, y + what.y, z + what.z);
    }

    double x;
    double y;
    double z;
};

}

namespace Racquetball{

class Ball{
public:
    Physics::Vector position;
    Physics::Vector move;
};

class Behavior{
public:
};

class Player{
public:
    Util::ReferenceCount<Behavior> behavior;
    Physics::Vector position;
    Physics::Vector move;
};

}

int main(){
    if (!al_init()){
        std::cout << "Could not initialize Allegro5" << std::endl;
        return -1;
    }
    al_set_new_display_flags(ALLEGRO_WINDOWED);

    ALLEGRO_DISPLAY * display = al_create_display(800, 600);
    al_destroy_display(display);
    return 0;
}
