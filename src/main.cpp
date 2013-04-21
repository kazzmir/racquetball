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

    Vector operator+(const Vector & what) const {
        return Vector(x + what.x, y + what.y, z + what.z);
    }

    Vector operator*(const double magnitude) const {
        return Vector(x * magnitude, y * magnitude, z * magnitude);
    }

    Vector operator-(const Vector & what) const {
        return Vector(x - what.x, y - what.y, z - what.z);
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

    void debug() const {
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
    lookTheta(180),
    lookPhi(0){
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

    const Physics::Vector & getPosition() const {
        return position;
    }

    /* theta and phi should be in degrees */
    Physics::Vector onSphere(double theta, double phi) const {
        double r_theta = theta * ALLEGRO_PI / 180.0;
        double r_phi = phi * ALLEGRO_PI / 180.0;
        double x = sin(r_theta) * cos(r_phi);
        double z = cos(r_theta) * cos(r_phi);
        double y = sin(r_phi);
        return Physics::Vector(x, y, z);
    }

    const Physics::Vector getLook() const {
        return onSphere(lookTheta, lookPhi);
    }

    const Physics::Vector getLookPerpendicular() const {
        return onSphere(lookTheta + 90, lookPhi);
    }

    /* Move the eye. x rotates theta, y rotates phi */
    void changeLook(double x, double y){
        double sensitivity = 6;
        lookTheta -= x / sensitivity;
        lookPhi -= y / sensitivity;
        if (lookPhi > 90){
            lookPhi = 90;
        }
        if (lookPhi < -90){
            lookPhi = -90;
        }

        // std::cout << "Theta: " << lookTheta << " Phi: " << lookPhi << std::endl;
    }

    void move(const Physics::Vector & much){
        this->position = this->position + much;
    }

    /* position in space */
    Physics::Vector position;
    /* velocity */
    Physics::Vector velocity;
    /* where the camera is looking */
    // Physics::Vector look;

    double lookTheta;
    double lookPhi;
};

void al_look_at_transform(ALLEGRO_TRANSFORM *transform, const Physics::Vector & look, const Physics::Vector & up){
    ALLEGRO_TRANSFORM tmp;

    // Physics::Vector f = (look - camera).normalize();
    Physics::Vector f = look.normalize();
    Physics::Vector s = f.cross(up);
    Physics::Vector u = s.cross(f);

    al_identity_transform(&tmp);

    tmp.m[0][0] = s.getX();
    tmp.m[0][1] = u.getX();
    tmp.m[0][2] = -f.getX();
    tmp.m[0][3] = 0;

    tmp.m[1][0] = s.getY();
    tmp.m[1][1] = u.getY();
    tmp.m[1][2] = -f.getY();
    tmp.m[1][3] = 0;

    tmp.m[2][0] = s.getZ();
    tmp.m[2][1] = u.getZ();
    tmp.m[2][2] = -f.getZ();
    tmp.m[2][3] = 0;

    tmp.m[3][0] = 0;
    tmp.m[3][1] = 0;
    tmp.m[3][2] = 0;
    tmp.m[3][3] = 1;

    al_compose_transform(transform, &tmp);
}

class Cube{
public:
    Cube(int size){
        points[0].x = -size / 2;
        points[0].y = -size / 2;
        points[0].z = -size / 2;

        points[1].x = size / 2;
        points[1].y = -size / 2;
        points[1].z = -size / 2;

        points[2].x = -size / 2;
        points[2].y = -size / 2;
        points[2].z = size / 2;

        points[3].x = size / 2;
        points[3].y = -size / 2;
        points[3].z = size / 2;

        points[4].x = -size / 2;
        points[4].y = size / 2;
        points[4].z = -size / 2;

        points[5].x = size / 2;
        points[5].y = size / 2;
        points[5].z = -size / 2;

        points[6].x = -size / 2;
        points[6].y = size / 2;
        points[6].z = size / 2;

        points[7].x = size / 2;
        points[7].y = size / 2;
        points[7].z = size / 2;

        for (int i = 0; i < 8; i++){
            points[i].u = 0;
            points[i].v = 0;
        }
    }

    void draw(double x, double y, double z, const ALLEGRO_COLOR & color, int * indicies){
        ALLEGRO_TRANSFORM transform;
        al_identity_transform(&transform);
        al_translate_transform_3d(&transform, x, y, z);
        const ALLEGRO_TRANSFORM * current = al_get_current_transform();
        ALLEGRO_TRANSFORM save;
        al_copy_transform(&save, current);
        al_compose_transform(&transform, current);
        al_use_transform(&transform);

        for (int i = 0; i < 4; i++){
            points[indicies[i]].color = color;
        }

        al_draw_indexed_prim(points, NULL, NULL, indicies, 4, ALLEGRO_PRIM_TRIANGLE_FAN);
        al_use_transform(&save);
    }

    void drawLeft(double x, double y, double z, const ALLEGRO_COLOR & color){
        int indicies[4] = {0, 2, 6, 4};
        draw(x, y, z, color, indicies);
    }

    void drawRight(double x, double y, double z, const ALLEGRO_COLOR & color){
        int indicies[4] = {1, 3, 7, 5};
        draw(x, y, z, color, indicies);
    }
    
    void drawTop(double x, double y, double z, const ALLEGRO_COLOR & color){
        int indicies[4] = {0, 1, 3, 2};
        draw(x, y, z, color, indicies);
    }
    
    void drawBottom(double x, double y, double z, const ALLEGRO_COLOR & color){
        int indicies[4] = {4, 5, 7, 6};
        draw(x, y, z, color, indicies);
    }
    
    void drawBack(double x, double y, double z, const ALLEGRO_COLOR & color){
        int indicies[4] = {0, 1, 5, 4};
        draw(x, y, z, color, indicies);
    }
    
    void drawFront(double x, double y, double z, const ALLEGRO_COLOR & color){
        int indicies[4] = {2, 3, 7, 6};
        draw(x, y, z, color, indicies);
    }

    void draw(double x, double y, double z, const ALLEGRO_COLOR & color){
        drawLeft(x, y, z, color);
        drawRight(x, y, z, color);
        drawTop(x, y, z, color);
        drawBottom(x, y, z, color);
        drawFront(x, y, z, color);
        drawBack(x, y, z, color);
    }

    ALLEGRO_VERTEX points[8];
};

void draw(const Camera & camera){
    al_clear_to_color(al_map_rgb_f(0, 0, 0));
    al_clear_depth_buffer(1);
    ALLEGRO_TRANSFORM transform;
    al_identity_transform(&transform);
    al_translate_transform_3d(&transform, -camera.positionX(), -camera.positionY(), -camera.positionZ());
    Physics::Vector up(0, 1, 0);
    // Physics::Vector cross = up.cross(camera.getLook()).normalize();
    // cross.debug();
    // std::cout << "radians " << acos(up.dot(camera.getLook())) << std::endl;
    // al_rotate_transform_3d(&transform, cross.getX(), cross.getY(), cross.getZ(), acos(up.dot(camera.getLook())));
    // al_rotate_transform_3d(&transform, 0, 1, 0, 0);
    al_look_at_transform(&transform, camera.getLook(), up);
    
    ALLEGRO_DISPLAY *display = al_get_current_display();
    int dw = al_get_display_width(display);
    int dh = al_get_display_height(display);
    // al_perspective_transform(&transform, -180 * dw / dh, -180, 180, 180 * dw / dh, 180, 3000);
    float distance = 10;

    float aspect = (float) dw / (float) dh;
    float fovx = 90;
    float fovy = fovx / aspect;
    /*
    float left = -dw / 4;
    float right = dw / 4;
    float top = -dh / 4;
    float bottom = dh / 4;
    */

    float left = -distance / 2 * aspect;
    float right = distance / 2 * aspect;
    float top = -distance / 2;
    float bottom = distance / 2;

    /*
    std::cout << "Fovy: " << fovy << std::endl;
    std::cout << "Y top: " << (1 / tan(fovy * ALLEGRO_PI / 180.0 / 2)) << std::endl;
    */

    // float right = tan(fov * ALLEGRO_PI / 180 / 2) * near * 2 + left;
    float near = (right - left) / 2 / tan(fovx * ALLEGRO_PI / 180 / 2);
    float far = near + 1000;

    /*
    std::cout << "Tan of fov: " << tan(fovx * ALLEGRO_PI / 180.0 / 2.0) << std::endl;
    std::cout << "Fov: " << fovx << std::endl;
    */
    // std::cout << "Near z " << near << std::endl;

    // al_perspective_transform(&transform, -distance * dw / dh, -distance, 50, distance * dw / dh, distance, 1000);
    al_perspective_transform(&transform, left, top, near, right, bottom, far);
    // al_perspective_transform(&transform, -1, -1, 5, 1, 1, far);

    al_set_projection_transform(display, &transform);

    al_set_blender(ALLEGRO_ADD, ALLEGRO_ALPHA, ALLEGRO_INVERSE_ALPHA);

    Cube cube(500);

    cube.drawLeft(0, 0, 0, al_map_rgba_f(1, 0.8, 0.3, 0.8));
    cube.drawRight(0, 0, 0, al_map_rgba_f(0.8, 1, 0.8, 0.8));
    cube.drawTop(0, 0, 0, al_map_rgba_f(1, 0.8, 1, 0.8));
    cube.drawBottom(0, 0, 0, al_map_rgba_f(0.3, 0.5, 0.8, 0.8));
    cube.drawBack(0, 0, 0, al_map_rgba_f(0.7, 0.4, 0.8, 0.8));
    cube.drawFront(0, 0, 0, al_map_rgba_f(0.2, 0.9, 0.5, 0.8));

    cube.draw(-600, 0, 0, al_map_rgb_f(1, 1, 1));
    cube.draw(600, 0, 0, al_map_rgb_f(1, 1, 1));
    
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
    al_install_mouse();
    al_init_primitives_addon();
    al_set_new_display_flags(ALLEGRO_WINDOWED);

    ALLEGRO_TIMER * timer = al_create_timer(ALLEGRO_BPS_TO_SECS(60));
    al_set_new_display_option(ALLEGRO_DEPTH_SIZE, 24, ALLEGRO_REQUIRE);
    ALLEGRO_DISPLAY * display = al_create_display(800, 600);
    al_set_render_state(ALLEGRO_DEPTH_TEST, 1);
    al_set_render_state(ALLEGRO_DEPTH_FUNCTION, ALLEGRO_RENDER_LESS);
    ALLEGRO_EVENT_QUEUE * queue = al_create_event_queue();
    al_register_event_source(queue, al_get_keyboard_event_source());
    al_register_event_source(queue, al_get_mouse_event_source());
    al_register_event_source(queue, al_get_timer_event_source(timer));

    al_hide_mouse_cursor(display);
    al_grab_mouse(display);
    Racquetball::Camera camera;

    camera.move(Physics::Vector(0, 0, 600));
    al_start_timer(timer);

    ALLEGRO_EVENT event;
    bool done = false;
    double speed = 1;
    while (!done){
        bool draw = false;
        while (al_get_next_event(queue, &event)){
            switch (event.type){
                case ALLEGRO_EVENT_TIMER: {
                    draw = true;
                    break;
                }
                case ALLEGRO_EVENT_MOUSE_AXES: {
                    camera.changeLook(event.mouse.dx, event.mouse.dy);
                    // camera.getLook().debug();
                    al_set_mouse_xy(display, al_get_display_width(display) / 2, al_get_display_height(display) / 2);
                    break;
                }
                case ALLEGRO_EVENT_KEY_CHAR: {
                    switch (event.keyboard.keycode){
                        case ALLEGRO_KEY_W:
                        case ALLEGRO_KEY_UP: {
                            // camera.move(Physics::Vector(0, 0, -speed));
                            camera.move(camera.getLook() * speed);
                            break;
                        }
                        case ALLEGRO_KEY_X:
                        case ALLEGRO_KEY_DOWN: {
                            // camera.move(Physics::Vector(0, 0, speed));
                            camera.move(camera.getLook() * -speed);
                            break;
                        }
                        case ALLEGRO_KEY_A:
                        case ALLEGRO_KEY_LEFT: {
                            camera.move(camera.getLookPerpendicular() * speed);
                            break;
                        }
                        case ALLEGRO_KEY_D:
                        case ALLEGRO_KEY_RIGHT: {
                            camera.move(camera.getLookPerpendicular() * -speed);
                            break;
                        }
                        /*
                        case ALLEGRO_KEY_Z: {
                            camera.move(Physics::Vector(0, speed, 0));
                            break;
                        }
                        case ALLEGRO_KEY_X: {
                            camera.move(Physics::Vector(0, -speed, 0));
                            break;
                        }
                        */
                        case ALLEGRO_KEY_ESCAPE: done = true; break;
                    }
                    camera.getPosition().debug();
                    break;
                }
            }

            // camera.getPosition().debug();
        }

        if (draw){
            Racquetball::draw(camera);
        } else {
            al_rest(0.001);
        }
    }

    /*
    for (int i = 0; i < 80; i++){
        Racquetball::draw(camera);
        camera.move(Physics::Vector(0, 0, -0.5));
        al_rest(0.1);
    }
    */

    al_ungrab_mouse();
    al_destroy_display(display);
    return 0;
}
