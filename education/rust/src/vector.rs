#[derive(Clone, Debug, Default)]
pub struct Vector<T> {
    data: Vec<T>, shape: Vec<usize>
}

impl<T: Clone + Default> Vector<T> {
    pub fn data(&mut self) -> &mut Vec<T> {
        return &mut self.data
    }
    pub fn init(shape: Vec<usize>) -> Vector<T> {
        return Vector {data: vec![T::default(); shape.iter().product()], shape: shape}
    }
    pub fn new(data: Vec<T>) -> Vector<T> {
        return Vector {shape: vec![data.len()], data: data}
    }
    pub fn shape(&self) -> &Vec<usize> {
        return &self.shape
    }
    pub fn size(&self) -> usize {
        return self.data.len()
    }
}
