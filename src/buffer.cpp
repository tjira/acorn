#include "buffer.h"

Buffer::~Buffer() {
    glDeleteVertexArrays(1, &vao), glDeleteBuffers(1, &vbo);
};

Buffer& Buffer::operator=(const Buffer& buffer) {
    glDeleteVertexArrays(1, &vao), glDeleteBuffers(1, &vbo); this->data = buffer.data, generate(); return *this;
}

void Buffer::bind() const {
    glBindVertexArray(vao);
}

void Buffer::generate() {
    // generate and bind the buffers
    glGenVertexArrays(1, &vao), glGenBuffers(1, &vbo), glBindBuffer(GL_ARRAY_BUFFER, vbo), glBindVertexArray(vao);

    // set position, normal and collor attributes
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, position));
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, normal));
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, color));

    // enable attributes
    glEnableVertexAttribArray(0), glEnableVertexAttribArray(1), glEnableVertexAttribArray(2);

    // set the data
    glBufferData(GL_ARRAY_BUFFER, data.size() * sizeof(Vertex), data.data(), GL_STATIC_DRAW);
}
