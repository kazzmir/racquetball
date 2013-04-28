#include <allegro5/allegro.h>
#include <allegro5/allegro_primitives.h>
#include <iostream>
#include <math.h>
#include <vector>
#include <time.h>
#include <stdlib.h>

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

    Vector operator*(const Vector & what) const {
        return Vector(x * what.x, y * what.y, z * what.z);
    }

    Vector operator/(const double magnitude) const {
        return *this * (1 / magnitude);
    }

    Vector operator-(const Vector & what) const {
        return Vector(x - what.x, y - what.y, z - what.z);
    }

    Vector operator-() const {
        return Vector(-x, -y, -z);
    }

    Vector & operator+=(const Vector & what){
        x += what.getX();
        y += what.getY();
        z += what.getZ();
        return *this;
    }

    Vector & operator-=(const Vector & what){
        x -= what.getX();
        y -= what.getY();
        z -= what.getZ();
        return *this;
    }
    
    void setZ(double z){
        this->z = z;
    }
    
    void setX(double x){
        this->x = x;
    }

    void setY(double y){
        this->y = y;
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

    double distance(const Vector & him) const {
        double xs = x - him.getX();
        double ys = y - him.getY();
        double zs = z - him.getZ();
        return sqrt(xs * xs + ys * ys + zs * zs);
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

Vector gravity(0, 0.7, 0);

}

namespace Racquetball{

class Translation{
public:
    Translation(double x, double y, double z){
        ALLEGRO_TRANSFORM transform;
        al_identity_transform(&transform);
        al_translate_transform_3d(&transform, x, y, z);
        const ALLEGRO_TRANSFORM * current = al_get_current_transform();
        al_copy_transform(&save, current);
        al_compose_transform(&transform, current);
        al_use_transform(&transform);
    }

    ~Translation(){
        al_use_transform(&save);
    }

    ALLEGRO_TRANSFORM save;
};

class Triangle{
public:
    Triangle(const Physics::Vector & point1, const Physics::Vector & point2, const Physics::Vector & point3):
    point1(point1),
    point2(point2),
    point3(point3){
    }

    std::vector<Triangle> subdivide(double radius) const {
        std::vector<Triangle> out;

        Physics::Vector mid12 = ((point1 + point2) / 2).normalize() * radius;
        Physics::Vector mid23 = ((point2 + point3) / 2).normalize() * radius;
        Physics::Vector mid13 = ((point1 + point3) / 2).normalize() * radius;

        out.push_back(Triangle(point1, mid12, mid13));
        out.push_back(Triangle(mid12, point2, mid23));
        out.push_back(Triangle(mid12, mid23, mid13));
        out.push_back(Triangle(mid13, mid23, point3));

        return out;
    }

    void setup_vertex(ALLEGRO_VERTEX * vertex, const Physics::Vector & point) const {
        vertex->x = point.getX();
        vertex->y = point.getY();
        vertex->z = point.getZ();
        vertex->u = 0;
        vertex->v = 0;
    }

    void draw(const ALLEGRO_COLOR & color) const {
        ALLEGRO_VERTEX points[3];
        setup_vertex(&points[0], point1);
        setup_vertex(&points[1], point2);
        setup_vertex(&points[2], point3);

        for (int i = 0; i < 3; i++){
            points[i].color = color;
        }

        al_draw_prim(points, NULL, NULL, 0, 4, ALLEGRO_PRIM_TRIANGLE_FAN);
        for (int i = 0; i < 3; i++){
            points[i].color = al_map_rgb_f(0, 0, 0);
        }
        al_draw_prim(points, NULL, NULL, 0, 4, ALLEGRO_PRIM_LINE_LOOP);
    }

    Physics::Vector point1, point2, point3;
};



class Sphere{
public:
    Sphere(double radius, int level){
        std::vector<Triangle> octohedron = initialOctohedron(radius);
        triangles = generate(octohedron, radius, level);
    }

    std::vector<Triangle> triangles;

    std::vector<Triangle> initialOctohedron(double radius){
        std::vector<Triangle> out;

        Physics::Vector xplus(1 * radius, 0, 0);
        Physics::Vector xminus(-1 * radius, 0, 0);
        Physics::Vector yplus(0, 1 * radius, 0);
        Physics::Vector yminus(0, -1 * radius, 0);
        Physics::Vector zplus(0, 0, 1 * radius);
        Physics::Vector zminus(0, 0, -1 * radius);
        
        out.push_back(Triangle(xplus, zplus, yplus));
        out.push_back(Triangle(yplus, zplus, xminus));
        out.push_back(Triangle(xminus, zplus, yminus));
        out.push_back(Triangle(yminus, zplus, xplus));
        out.push_back(Triangle(xplus, yplus, zminus));
        out.push_back(Triangle(yplus, xminus, zminus));
        out.push_back(Triangle(xminus, yminus, zminus));
        out.push_back(Triangle(yminus, xplus, zminus));

        return out;
    }

    std::vector<Triangle> generate(const std::vector<Triangle> & in, double radius, int level){
        if (level == 0){
            return in;
        }

        std::vector<Triangle> out;

        for (std::vector<Triangle>::const_iterator it = in.begin(); it != in.end(); it++){
            const Triangle & original = *it;
            std::vector<Triangle> more = original.subdivide(radius);
            out.insert(out.end(), more.begin(), more.end());
        }

        return generate(out, radius, level - 1);
    }

    void draw(double x, double y, double z, const ALLEGRO_COLOR & color) const {
        Translation translate(x, y, z);

        for (std::vector<Triangle>::const_iterator it = triangles.begin(); it != triangles.end(); it++){
            const Triangle & triangle = *it;
            triangle.draw(color);
        }
    }
};

class Ball{
public:
    Ball(const Physics::Vector & position):
    position(position),
    velocity(0, 0, 0),
    model(10, 3){
    }

    double getSize() const {
        return 10;
    }

    void setPosition(const Physics::Vector & what){
        this->position = what;
    }

    void move(){
        position += velocity;
        velocity -= Physics::gravity;
    }

    const Physics::Vector & getPosition() const {
        return position;
    }

    void setVelocity(const Physics::Vector & velocity){
        this->velocity = velocity;
    }

    const Physics::Vector & getVelocity() const {
        return velocity;
    }

    void draw(double x, double y, double z) const {
        Translation translate(x, y, z);
        model.draw(position.getX(), position.getY(), position.getZ(), al_map_rgb_f(0.3, 0.2, 1));
    }

    /* Sort of models different ball types.
     * pink = fast
     * green = normal
     * black = slow
     */
    double density;
    Physics::Vector position;
    Physics::Vector velocity;
    Sphere model;
};

class Cube{
public:
    Cube(int size):
    size(size){
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

    int size;

    double getWidth() const {
        return size;
    }

    double getHeight() const {
        return size;
    }

    double getLength() const {
        return size;
    }

    void draw(double x, double y, double z, const ALLEGRO_COLOR & color, int * indicies) const {
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
        for (int i = 0; i < 4; i++){
            points[indicies[i]].color = al_map_rgb_f(0, 0, 0);
        }
        al_draw_indexed_prim(points, NULL, NULL, indicies, 4, ALLEGRO_PRIM_LINE_LOOP);
        al_use_transform(&save);
    }

    void drawLeft(double x, double y, double z, const ALLEGRO_COLOR & color) const {
        int indicies[4] = {0, 2, 6, 4};
        draw(x, y, z, color, indicies);
    }

    void drawRight(double x, double y, double z, const ALLEGRO_COLOR & color) const {
        int indicies[4] = {1, 3, 7, 5};
        draw(x, y, z, color, indicies);
    }
    
    void drawTop(double x, double y, double z, const ALLEGRO_COLOR & color) const {
        int indicies[4] = {4, 5, 7, 6};
        draw(x, y, z, color, indicies);
    }
    
    void drawBottom(double x, double y, double z, const ALLEGRO_COLOR & color) const {
        int indicies[4] = {0, 1, 3, 2};
        draw(x, y, z, color, indicies);
    }
    
    void drawBack(double x, double y, double z, const ALLEGRO_COLOR & color) const {
        int indicies[4] = {2, 3, 7, 6};
        draw(x, y, z, color, indicies);
    }
    
    void drawFront(double x, double y, double z, const ALLEGRO_COLOR & color) const {
        int indicies[4] = {0, 1, 5, 4};
        draw(x, y, z, color, indicies);
    }

    void draw(double x, double y, double z, const ALLEGRO_COLOR & color) const {
        drawLeft(x, y, z, color);
        drawRight(x, y, z, color);
        drawTop(x, y, z, color);
        drawBottom(x, y, z, color);
        drawFront(x, y, z, color);
        drawBack(x, y, z, color);
    }

    mutable ALLEGRO_VERTEX points[8];
};

class Behavior{
public:
};

class Player{
public:
    Player(const Physics::Vector & position):
    speed(15),
    position(position),
    move(0, 0, 0),
    model(10){
    }

    void moveLeft(const Physics::Vector & look){
        position += Physics::Vector(-look.getZ(), 0, look.getX()) * -speed;
    }

    void moveRight(const Physics::Vector & look){
        position += Physics::Vector(-look.getZ(), 0, look.getX()) * speed;
    }

    void moveForward(const Physics::Vector & look){
        position += Physics::Vector(look.getX(), 0, look.getZ()) * speed;
    }

    void moveBackward(const Physics::Vector & look){
        position += Physics::Vector(look.getX(), 0, look.getZ()) * -speed;
    }

    const Physics::Vector & getPosition() const {
        return position;
    }

    void draw(double x, double y, double z) const {
        Translation translate(x, y, z);
        model.draw(position.getX(), position.getY(), position.getZ(), al_map_rgb_f(1, 1, 1));
    }

    /* How fast the player can move */
    double speed;

    /* How hard the player can hit */
    double power;

    /* How accurate the player is */
    double talent;

    Util::ReferenceCount<Behavior> behavior;
    Physics::Vector position;
    Physics::Vector move;
    Cube model;
};

class Camera{
public:
    Camera(){
    }

    virtual ~Camera(){
    }

    virtual double positionX() const = 0;
    virtual double positionY() const = 0;
    virtual double positionZ() const = 0;

    virtual const Physics::Vector & getPosition() const = 0;
    virtual const Physics::Vector getLook() const = 0;
};

class StationaryCamera: public Camera {
public:
    StationaryCamera(const Physics::Vector & position, const Physics::Vector & look):
    position(position),
    look(look){
    }

    const Physics::Vector position;
    const Physics::Vector look;

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

    const Physics::Vector getLook() const {
        return look;
    }
};

class FreeCamera: public Camera {
public:
    FreeCamera():
    position(0, 0, 0),
    velocity(0, 0, 0),
    lookTheta(180),
    lookPhi(0){
    }

    /* position in space */
    Physics::Vector position;
    /* velocity */
    Physics::Vector velocity;
    /* where the camera is looking */
    // Physics::Vector look;

    double lookTheta;
    double lookPhi;

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

    /* where should be normalized */
    void lookAt(const Physics::Vector & where){
        double r_phi = asin(where.getY());
        double r_theta = asin(where.getX() / cos(r_phi));
        /* asin only returns a value in -pi/2,pi/2 so we have to move it
         * to quadrant 2/3
         */
        if (where.getZ() < 0){
            double r180 = ALLEGRO_PI;
            r_theta = r180 - r_theta;
        }
        lookTheta = r_theta * 180.0 / ALLEGRO_PI;
        lookPhi = r_phi * 180.0 / ALLEGRO_PI;
        std::cout << "Look at ";
        where.debug();
        std::cout << "Theta " << lookTheta << " Phi " << lookPhi << std::endl;
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
        if (lookTheta < 0){
            lookTheta += 360;
        }
        if (lookTheta > 360){
            lookTheta -= 360;
        }
        lookPhi += y / sensitivity;
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
};

void al_look_at_transform(ALLEGRO_TRANSFORM *transform, const Physics::Vector & look, const Physics::Vector & up){
    ALLEGRO_TRANSFORM tmp;

    // Physics::Vector f = (look - camera).normalize();
    Physics::Vector front = look.normalize();
    Physics::Vector side = front.cross(up).normalize();
    Physics::Vector cameraUp = side.cross(front);

    al_identity_transform(&tmp);

    tmp.m[0][0] = side.getX();
    tmp.m[0][1] = cameraUp.getX();
    tmp.m[0][2] = -front.getX();
    tmp.m[0][3] = 0;

    tmp.m[1][0] = side.getY();
    tmp.m[1][1] = cameraUp.getY();
    tmp.m[1][2] = -front.getY();
    tmp.m[1][3] = 0;

    tmp.m[2][0] = side.getZ();
    tmp.m[2][1] = cameraUp.getZ();
    tmp.m[2][2] = -front.getZ();
    tmp.m[2][3] = 0;

    tmp.m[3][0] = 0;
    tmp.m[3][1] = 0;
    tmp.m[3][2] = 0;
    tmp.m[3][3] = 1;

    al_compose_transform(transform, &tmp);
}

class Court{
public:
    Court():
    court(2000),
    ball(Physics::Vector(0, 50, 0)),
    player(Physics::Vector(-100, -court.getHeight() / 2 + 250, -100)){
        ball.setVelocity(Physics::Vector(0, 0, 0));
    }

    Cube court;

    double height;
    double length;
    double width;

    Ball ball;
    Player player;

    Player & getPlayer(){
        return player;
    }

    const Player & getPlayer() const {
        return player;
    }

    const Ball & getBall() const {
        return ball;
    }

    bool outOfBounds(const Physics::Vector & position){
        return position.getX() < -court.getWidth() / 2 ||
               position.getX() > court.getWidth() / 2 ||
               position.getY() < -court.getHeight() / 2 ||
               position.getY() > court.getHeight() / 2 ||
               position.getZ() < -court.getLength() / 2 ||
               position.getZ() > court.getLength() / 2;
    }

    void logic(){
        ball.move();

        double friction = 0.8;

        if (ball.getPosition().getX() - ball.getSize() < -court.getWidth() / 2){
            Physics::Vector where = ball.getPosition();
            where.setX(-court.getWidth() / 2 + ball.getSize());
            ball.setPosition(where);
            ball.setVelocity(ball.getVelocity() * Physics::Vector(-1, 1, 1) * friction);
        }

        if (ball.getPosition().getX() + ball.getSize() > court.getWidth() / 2){
            Physics::Vector where = ball.getPosition();
            where.setX(court.getWidth() / 2 - ball.getSize());
            ball.setPosition(where);
            ball.setVelocity(ball.getVelocity() * Physics::Vector(-1, 1, 1) * friction);
        }

        if (ball.getPosition().getZ() - ball.getSize() < -court.getLength() / 2){
            Physics::Vector where = ball.getPosition();
            where.setZ(-court.getLength() / 2 + ball.getSize());
            ball.setPosition(where);
            ball.setVelocity(ball.getVelocity() * Physics::Vector(1, 1, -1) * friction);
        }

        if (ball.getPosition().getZ() + ball.getSize() > court.getLength() / 2){
            Physics::Vector where = ball.getPosition();
            where.setZ(court.getLength() / 2 - ball.getSize());
            ball.setPosition(where);
            ball.setVelocity(ball.getVelocity() * Physics::Vector(1, 1, -1) * friction);
        }

        if (ball.getPosition().getY() - ball.getSize() < -court.getHeight() / 2){
            Physics::Vector where = ball.getPosition();
            where.setY(-court.getHeight() / 2 + ball.getSize());
            ball.setPosition(where);
            ball.setVelocity(ball.getVelocity() * Physics::Vector(1, -1, 1) * friction);
        }

        if (ball.getPosition().getY() + ball.getSize() > court.getHeight() / 2){
            Physics::Vector where = ball.getPosition();
            where.setY(court.getHeight() / 2 - ball.getSize());
            ball.setPosition(where);
            ball.setVelocity(ball.getVelocity() * Physics::Vector(1, -1, 1) * friction);
        }

        /*
        if (outOfBounds(ball.getPosition())){
            Physics::Vector next = -ball.getVelocity();
            if (ball.getPosition().getY() < -court.getHeight() / 2){
                Physics::Vector where = ball.getPosition();
                where.setY(-court.getHeight() / 2 + ball.getSize());
                ball.setPosition(where);
                next = next / 2;
            }
            ball.setVelocity(next);
        }
        */
    }

    double randomFloat(double what){
        return (double) (rand() % 1000) * what / 1000.0;
    }

    double randomFloat(double low, double high){
        return randomFloat(high - low) + low;
    }

    void hit(){
        ball.setVelocity(Physics::Vector(randomFloat(-5, 5) * 3, randomFloat(10) + 10, -40));
    }

    void draw(double x, double y, double z) const {
        court.drawLeft(x, y, z, al_map_rgb_f(0.5, 0.5, 0.5));
        court.drawRight(x, y, z, al_map_rgb_f(0.5, 0.5, 0.5));
        court.drawTop(x, y, z, al_map_rgb_f(1, 1, 0));
        court.drawFront(x, y, z, al_map_rgb_f(1, 1, 1));
        court.drawBack(x, y, z, al_map_rgb_f(1, 0, 0));
        court.drawBottom(x, y, z, al_map_rgb_f(0.8, 0.8, 0));
        ball.draw(x, y, z);
        player.draw(x, y, z);
    }

    /* Center court position (y = 0) */
    Physics::Vector getCenter(){
        return Physics::Vector(0, 0, 0);
    }
};

static void setup_draw_transform(const Camera & camera){
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
    float top = distance / 2;
    float bottom = -distance / 2;

    /*
    std::cout << "Fovy: " << fovy << std::endl;
    std::cout << "Y top: " << (1 / tan(fovy * ALLEGRO_PI / 180.0 / 2)) << std::endl;
    */

    // float right = tan(fov * ALLEGRO_PI / 180 / 2) * near * 2 + left;
    float near = (right - left) / 2 / tan(fovx * ALLEGRO_PI / 180 / 2);
    float far = near + 3000;

    /*
    std::cout << "Tan of fov: " << tan(fovx * ALLEGRO_PI / 180.0 / 2.0) << std::endl;
    std::cout << "Fov: " << fovx << std::endl;
    */
    // std::cout << "Near z " << near << std::endl;

    // al_perspective_transform(&transform, -distance * dw / dh, -distance, 50, distance * dw / dh, distance, 1000);
    al_perspective_transform(&transform, left, top, near, right, bottom, far);
    // al_perspective_transform(&transform, -1, -1, 5, 1, 1, far);

    al_set_projection_transform(display, &transform);
}

StationaryCamera computeStationaryCamera(const Physics::Vector & cameraLook, const Ball & ball, const Player & player){
    /*
    Physics::Vector direction = (player.getPosition() - ball.getPosition()).normalize();

    Physics::Vector position = player.getPosition() + direction * 15 + Physics::Vector(0, 10, 0);
    Physics::Vector look = (-(position - ball.getPosition())).normalize();
    // Physics::Vector look = (player.getPosition() - ball.getPosition()).normalize();
    */
    Physics::Vector position = player.getPosition() - cameraLook * 22 + Physics::Vector(0, 10, 0);
    return StationaryCamera(position, cameraLook);
}

void draw(const Court & court, const Camera & camera){
    al_clear_to_color(al_map_rgb_f(0, 0, 0));
    al_clear_depth_buffer(1);

    StationaryCamera stationary = computeStationaryCamera(camera.getLook(), court.getBall(), court.getPlayer());

    // setup_draw_transform(camera);
    setup_draw_transform(stationary);

    al_set_blender(ALLEGRO_ADD, ALLEGRO_ALPHA, ALLEGRO_INVERSE_ALPHA);

    court.draw(0, 0, 0);

    /*
    Cube cube(1000);

    cube.drawLeft(0, 0, 0, al_map_rgb_f(1, 0.8, 0.3));
    cube.drawRight(0, 0, 0, al_map_rgb_f(0.8, 1, 0.8));
    cube.drawTop(0, 0, 0, al_map_rgb_f(1, 0.8, 1));
    cube.drawBottom(0, 0, 0, al_map_rgb_f(0.3, 0.5, 0.8));
    cube.drawBack(0, 0, 0, al_map_rgb_f(0.7, 0.4, 0.8));
    cube.drawFront(0, 0, 0, al_map_rgb_f(0.2, 0.9, 0.5));
    */

    /*
    cube.draw(-600, 0, 0, al_map_rgb_f(1, 1, 1));
    cube.draw(600, 0, 0, al_map_rgb_f(1, 1, 1));
    */
    
    /*
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
    */
    /*
    Sphere sphere(40, 0);
    sphere.draw(0, 0, 0, al_map_rgb_f(1, 0, 0));
    
    Sphere sphere1(40, 1);
    sphere1.draw(120, 0, 0, al_map_rgb_f(0, 1, 0));

    Sphere sphere2(40, 2);
    sphere2.draw(210, 0, 0, al_map_rgb_f(0, 0, 1));
    
    Sphere sphere3(40, 3);
    sphere3.draw(300, 0, 0, al_map_rgb_f(1, 1, 0));
    
    Sphere sphere4(40, 4);
    sphere4.draw(400, 0, 0, al_map_rgb_f(1, 1, 1));
    */

    al_flip_display();
}

}

static bool init_allegro(){
    if (ALLEGRO_VERSION_INT != al_get_allegro_version()){
        std::cout << "Wrong Allegro5 version. Compiled with " << ALLEGRO_VERSION_INT << " but running with " << al_get_allegro_version() << std::endl;
        return false;
    }
    if (!al_init()){
        std::cout << "Could not initialize Allegro5" << std::endl;
        return false;
    }
    if (!al_install_keyboard()){
        std::cout << "Could not initialize keyboard" << std::endl;
        return false;
    }
    if (!al_install_mouse()){
        std::cout << "Could not initialize mouse" << std::endl;
        return false;
    }
    if (!al_init_primitives_addon()){
        std::cout << "Could not initialize primitives addon" << std::endl;
        return false;
    }

    return true;
}

ALLEGRO_DISPLAY * setup_display(int width, int height){
    al_set_new_display_flags(ALLEGRO_WINDOWED);
    al_set_new_display_option(ALLEGRO_DEPTH_SIZE, 24, ALLEGRO_REQUIRE);
    ALLEGRO_DISPLAY * display = al_create_display(width, height);
    /* Enable the depth buffer */
    al_set_render_state(ALLEGRO_DEPTH_TEST, 1);
    al_set_render_state(ALLEGRO_DEPTH_FUNCTION, ALLEGRO_RENDER_LESS);
    return display;
}

int main(){
    if (!init_allegro()){
        return 1;
    }

    srand(time(NULL));

    ALLEGRO_DISPLAY * display = setup_display(1024, 768);

    ALLEGRO_EVENT_QUEUE * queue = al_create_event_queue();

    ALLEGRO_TIMER * timer = al_create_timer(ALLEGRO_BPS_TO_SECS(60));
    al_register_event_source(queue, al_get_keyboard_event_source());
    al_register_event_source(queue, al_get_mouse_event_source());
    al_register_event_source(queue, al_get_timer_event_source(timer));

    al_hide_mouse_cursor(display);
    al_grab_mouse(display);
    Racquetball::FreeCamera camera;

    camera.move(Physics::Vector(0, 0, 200));
    al_start_timer(timer);

    Racquetball::Court court;

    struct Hold{
        Hold():
        left(false),
        right(false),
        up(false),
        down(false),
        left_click(false),
        right_click(false){
        }

        bool left, right, up, down;
        bool left_click, right_click;
    };

    Hold hold;
    bool hit = false;

    ALLEGRO_EVENT event;
    bool done = false;
    double speed = 3;
    while (!done){
        bool draw = false;
        while (al_get_next_event(queue, &event)){
            switch (event.type){
                case ALLEGRO_EVENT_TIMER: {
                    Racquetball::Player & player = court.getPlayer();
                    if (hold.up){
                        player.moveForward(camera.getLook());
                        camera.move(camera.getLook() * speed);
                    }
                    if (hold.down){
                        player.moveBackward(camera.getLook());
                        camera.move(camera.getLook() * -speed);
                    }
                    if (hold.left){
                        player.moveLeft(camera.getLook());
                        camera.move(camera.getLookPerpendicular() * speed);
                    }
                    if (hold.right){
                        player.moveRight(camera.getLook());
                        camera.move(camera.getLookPerpendicular() * -speed);
                    }

                    if (hold.right_click){
                        camera.lookAt((court.getBall().getPosition() - court.getPlayer().getPosition()).normalize());
                    }

                    if (hold.left_click && court.getPlayer().getPosition().distance(court.getBall().getPosition()) < 300){
                        court.hit();
                        hit = false;
                    }

                    hold.left_click = false;

                    court.logic();

                    draw = true;
                    break;
                }
                case ALLEGRO_EVENT_MOUSE_BUTTON_UP:
                case ALLEGRO_EVENT_MOUSE_BUTTON_DOWN: {
                    bool pressed = event.type == ALLEGRO_EVENT_MOUSE_BUTTON_DOWN;
                    if (event.mouse.button == 1){
                        hold.left_click = pressed;
                    }
                    if (event.mouse.button == 2){
                        hold.right_click = pressed;
                    }
                    break;
                }
                case ALLEGRO_EVENT_MOUSE_AXES: {
                    camera.changeLook(event.mouse.dx, event.mouse.dy);
                    // camera.getLook().debug();
                    /* Recenter mouse so we can get more delta mouse movements */
                    al_set_mouse_xy(display, al_get_display_width(display) / 2, al_get_display_height(display) / 2);
                    break;
                }
                case ALLEGRO_EVENT_KEY_DOWN:
                case ALLEGRO_EVENT_KEY_UP: {
                    bool pressed = event.type == ALLEGRO_EVENT_KEY_DOWN;

                    switch (event.keyboard.keycode){
                        case ALLEGRO_KEY_SPACE: {
                            if (pressed){
                                hit = true;
                            }
                            break;
                        }
                        case ALLEGRO_KEY_W:
                        case ALLEGRO_KEY_UP: {
                            hold.up = pressed;
                            break;
                        }
                        case ALLEGRO_KEY_X:
                        case ALLEGRO_KEY_DOWN: {
                            hold.down = pressed;
                            break;
                        }
                        case ALLEGRO_KEY_A:
                        case ALLEGRO_KEY_LEFT: {
                            hold.left = pressed;
                            break;
                        }
                        case ALLEGRO_KEY_D:
                        case ALLEGRO_KEY_RIGHT: {
                            hold.right = pressed;
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
                    break;
                }
            }

            // camera.getPosition().debug();
        }

        if (draw){
            Racquetball::draw(court, camera);
        } else {
            al_rest(0.001);
        }
    }

    al_ungrab_mouse();
    al_destroy_display(display);
    return 0;
}
