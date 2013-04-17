#include <allegro5/allegro.h>
#include <iostream>

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
