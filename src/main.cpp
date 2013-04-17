#include <allegro5/allegro.h>
#include <allegro5/allegro_primitives.h>
#include <iostream>
#include <math.h>

#include "pointer.h"

namespace Physics{

/* 3d vector.
 * x = left/right
 * y = up/down
 * z = forward/back
 */
class Vector{
public:
    Vector(double x, double y, double z):
    x(x), y(y), z(z){
    }

    Vector operator+(const Vector & what){
        return Vector(x + what.x, y + what.y, z + what.z);
    }

    double getX() const {
        return x;
    }
    
    double getY() const {
        return y;
    }
    
    double getZ() const {
        return z;
    }

    void debug(){
        std::cout << "x: " << x << " y: " << y << " z: " << z << std::endl;
    }

    Vector normalize() const {
        double large = x * x + y * y + z * z;
        if (large == 0){
            return Vector(0, 0, 0);
        }
        double magnitude = sqrt(large);
        return Vector(x / magnitude, y / magnitude, z / magnitude);
    }

    double dot(const Vector & what) const {
        return x * what.x + y * what.y + z * what.z;
    }

    Vector cross(const Vector & what){
        return Vector(y * what.z - z * what.y,
                      z * what.x - x * what.z,
                      x * what.y - y * what.x);
    }

protected:
    double x;
    double y;
    double z;
};

}

namespace Racquetball{

class Ball{
public:
    /* Sort of models different ball types.
     * pink = fast
     * green = normal
     * black = slow
     */
    double density;
    Physics::Vector position;
    Physics::Vector move;
};

class Behavior{
public:
};

class Player{
public:
    /* How fast the player can move */
    double speed;

    /* How hard the player can hit */
    double power;

    /* How accurate the player is */
    double talent;

    Util::ReferenceCount<Behavior> behavior;
    Physics::Vector position;
    Physics::Vector move;
};

class Court{
public:
    /* Center court position (y = 0) */
    Physics::Vector getCenter();

    double height;
    double length;
    double width;
};

class Camera{
public:
    Camera():
    position(0, 0, 0),
    velocity(0, 0, 0),
    look(0, 0, -1){
    }

    double positionX() const {
        return position.getX();
    }
    
    double positionY() const {
        return position.getY();
    }
    
    double positionZ() const {
        return position.getZ();
    }

    const Physics::Vector & getLook() const {
        return look;
    }

    void move(const Physics::Vector & much){
        this->position = this->position + much;
    }

    /* position in space */
    Physics::Vector position;
    /* velocity */
    Physics::Vector velocity;
    /* where the camera is looking */
    Physics::Vector look;
};

void draw(const Camera & camera){
    al_clear_to_color(al_map_rgb_f(0, 0, 0));
    ALLEGRO_TRANSFORM transform;
    al_identity_transform(&transform);
    al_translate_transform_3d(&transform, -camera.positionX(), -camera.positionY(), -camera.positionZ());
    Physics::Vector up(0, 1, 0);
    Physics::Vector cross = up.cross(camera.getLook()).normalize();
    cross.debug();
    std::cout << "radians " << acos(up.dot(camera.getLook())) << std::endl;
    // al_rotate_transform_3d(&transform, cross.getX(), cross.getY(), cross.getZ(), acos(up.dot(camera.getLook())));
    al_rotate_transform_3d(&transform, 0, 1, 0, 0);
    
    ALLEGRO_DISPLAY *display = al_get_current_display();
    int dw = al_get_display_width(display);
    int dh = al_get_display_height(display);
    // al_perspective_transform(&transform, -180 * dw / dh, -180, 180, 180 * dw / dh, 180, 3000);
    int distance = 10;
    al_perspective_transform(&transform, -distance * dw / dh, -distance, 50, distance * dw / dh, distance, 1000);
    al_set_projection_transform(display, &transform);

    /*
    for (int i = -100; i < 100; i++){
        ALLEGRO_VERTEX v[1];
        v[0].x = i;
        v[0].y = i;
        v[0].z = -10;
        v[0].color = al_map_rgb_f(1, 1, 1);
        v[0].u = 128;
        v[0].v = 0;
        al_draw_prim(v, NULL, NULL, 0, 1, ALLEGRO_PRIM_POINT_LIST);
    }
    */

    /*
    for (int y = 0; y < 20; y++){
        for (int x = 0; x < 360; x += 5){
            ALLEGRO_VERTEX v[1];
            v[0].x = cos(x * ALLEGRO_PI / 180.0) * (20 - y);
            v[0].y = y * 2;
            v[0].z = -4 + sin(x * ALLEGRO_PI / 180) * (20 - y);
            v[0].color = al_map_rgb_f(1, 1, 1);
            v[0].u = 128;
            v[0].v = 0;
            al_draw_prim(v, NULL, NULL, 0, 1, ALLEGRO_PRIM_POINT_LIST);

            v[0].x = cos(x * ALLEGRO_PI / 180.0) * (20 - y);
            v[0].y = -y * 2;
            v[0].z = -4 + sin(x * ALLEGRO_PI / 180) * (20 - y);
            v[0].color = al_map_rgb_f(1, 1, 1);
            v[0].u = 128;
            v[0].v = 0;
            al_draw_prim(v, NULL, NULL, 0, 1, ALLEGRO_PRIM_POINT_LIST);

        }
    }
    */

    double radius = 20;
    for (int phi = 0; phi < 360; phi += 5){
        for (int rho = 0; rho < 360; rho += 5){
            ALLEGRO_VERTEX v[1];
            double r_phi = phi * ALLEGRO_PI / 180.0;
            double r_rho = rho * ALLEGRO_PI / 180.0;
            v[0].x = radius * sin(r_phi) * cos(r_rho);
            v[0].y = radius * sin(r_phi) * sin(r_rho);
            v[0].z = radius * cos(r_phi);
            v[0].color = al_map_rgb_f(fabs(cos(r_phi)), fabs(cos(r_rho)), 1);
            v[0].u = 0;
            v[0].v = 0;
            al_draw_prim(v, NULL, NULL, 0, 1, ALLEGRO_PRIM_POINT_LIST);
        }
    }

    al_flip_display();
}

}

int main(){
    if (ALLEGRO_VERSION_INT != al_get_allegro_version()){
        std::cout << "Wrong Allegro5 version. Compiled with " << ALLEGRO_VERSION_INT << " but running with " << al_get_allegro_version() << std::endl;
    }
    if (!al_init()){
        std::cout << "Could not initialize Allegro5" << std::endl;
        return -1;
    }
    al_install_keyboard();
    al_init_primitives_addon();
    al_set_new_display_flags(ALLEGRO_WINDOWED);

    ALLEGRO_DISPLAY * display = al_create_display(800, 600);
    ALLEGRO_EVENT_QUEUE * queue = al_create_event_queue();
    al_register_event_source(queue, al_get_keyboard_event_source());

    Racquetball::Camera camera;

    camera.move(Physics::Vector(0, 0, 120));

    ALLEGRO_EVENT event;
    bool done = false;
    bool first = true;
    double speed = 0.5;
    while (!done){
        bool display = false;
        while (al_get_next_event(queue, &event)){
            switch (event.type){
                case ALLEGRO_EVENT_KEY_CHAR: {
                    std::cout << "Pressed " << event.keyboard.keycode << std::endl;
                    switch (event.keyboard.keycode){
                        case 84: {
                            display = true;
                            camera.move(Physics::Vector(0, 0, -speed));
                            break;
                        }
                        case 85: {
                            display = true;
                            camera.move(Physics::Vector(0, 0, speed));
                            break;
                        }
                        case 59: done = true; break;
                    }
                    break;
                }
            }
        }

        if (first || display){
            Racquetball::draw(camera);
        } else {
            al_rest(0.01);
        }

        first = false;
    }

    /*
    for (int i = 0; i < 80; i++){
        Racquetball::draw(camera);
        camera.move(Physics::Vector(0, 0, -0.5));
        al_rest(0.1);
    }
    */

    al_destroy_display(display);
    return 0;
}
