#include "geometry.h"

Geometry Geometry::Load(std::stringstream& file) {
    // define the molecule, length and containers for strings and extract the header
    Geometry molecule; int length;  std::string line, name; std::getline(file, line), std::getline(file, name), length = std::stoi(line);

    // extract first coords line
    std::getline(file, line);

    // add atom for each line.
    for (int i = 0; i < length; std::getline(file, line), i++) {

        // extract the coordinates
        std::string atom; float x, y, z; std::stringstream iss(line); iss >> atom >> x >> y >> z;

        // define the scale matrix and translate matrix
        glm::mat4 scale = glm::scale(glm::mat4(1), glm::vec3(ATOMSIZEFACTOR * ptable.at(atom).radius)), translate = glm::translate(glm::mat4(1.0f), glm::vec3(x, y, z));

        // add the atom to the molecule
        molecule.objects.push_back({translate, glm::mat4(1.0f), scale, atom});
     }

    // add bonds and return
    molecule.rebind(BINDINGFACTOR); return molecule;
}

glm::vec3 Geometry::get_center() const {
    // define the center and number of atoms
    glm::vec3 center(0); float size = 0;

    // sum all the coordinates
    for (const Object& object : objects) if (object.name != "bond" && object.name != "El") center += object.get_position(), size += 1;

    // return the center
    return center / size;
}

void Geometry::move(const glm::vec3& vector) {
    for (Object& object : objects) object.translate = glm::translate(object.translate, vector);
}

void Geometry::rebind(float factor) {
    // define the bond size and the new array of objects
    static float bond_size = BONDSIZE; std::vector<Object> objects; int atoms = 0;

    // add all the atoms to the new vector and extract the bond sizes from the previous bonds
    for (Object obj : this->objects) obj.name == "bond" ? (void)(bond_size = obj.scale[0][0]) : objects.push_back(obj), atoms++;

    // loop over all pairs of atoms
    for (int i = 0; i < atoms; i++) for (int j = i + 1; j < atoms; j++) {

        // extract the distances between the two atoms
        float distance = glm::length(objects.at(j).get_position() - objects.at(i).get_position());

        // ignore electrons
        if (objects.at(i).name == "El" || objects.at(j).name == "El") continue;

        // if the bond should be there
        if (distance < factor * (ptable.at(objects.at(i).name).covalent + ptable.at(objects.at(j).name).covalent)) {

            // extract the bond position and the bond vector
            glm::vec3 position = (objects.at(i).get_position() + objects.at(j).get_position()) / 2.0f, vector = objects.at(j).get_position() - objects.at(i).get_position();

            // calculate the cross product with the up vector and the angle of the bond with the up vctor
            glm::vec3 cross = glm::cross(glm::vec3(0, 1, 0), vector); float angle = atan2f(glm::length(cross), glm::dot(glm::vec3(0, 1, 0), vector));

            // scale, rotate and translate the bond
            glm::mat4 scale = glm::scale(glm::mat4(1), {bond_size, glm::length(vector) / 2.0f, bond_size}), rotate = glm::rotate(glm::mat4(1), angle, glm::normalize(cross)), translate = glm::translate(glm::mat4(1.0f), position);

            // add the bond to the objects
            objects.push_back({translate, rotate, scale, "bond"});
        }
    }

    // change the objects
    this->objects = objects;
};

void Geometry::render(const Shader& shader) const {
    // loop over every object
    for (size_t i = 0; i < objects.size(); i++) {

        // render bonds and atoms
        if (objects.at(i).name == "bond") meshes.at("bond")            .render(shader, objects.at(i).get_model());
        else                              meshes.at(objects.at(i).name).render(shader, objects.at(i).get_model());
    }
}
