#include "expression.h"

Expression::Expression(const std::string& expression_string, const std::vector<std::string>& variable_strings) : variables(variable_strings.size()) {
    // set the expression string and variable strings
    this->expression_string = expression_string; this->variable_strings = variable_strings;

    // define the symbol table
    exprtk::symbol_table<double> symbols;

    // add variables and constants to the expression
    for (size_t i = 0; i < variable_strings.size(); i++) {
        symbols.add_variable(variable_strings.at(i), variables.at(i));
    } symbols.add_constants();

    // register the symbol table
    expression.register_symbol_table(symbols);

    // parse the expression
    exprtk::parser<double>().compile(expression_string, expression);
}

Expression::Expression(const Expression& other) : variables(other.variables.size()) {
    // set the expression string and variable strings
    this->expression_string = other.expression_string; this->variable_strings = other.variable_strings;

    // define the symbol table
    exprtk::symbol_table<double> symbols;

    // add variables and constants to the expression
    for (size_t i = 0; i < variable_strings.size(); i++) {
        symbols.add_variable(variable_strings.at(i), variables.at(i));
    } symbols.add_constants();

    // register the symbol table
    expression.register_symbol_table(symbols);

    // parse the expression
    exprtk::parser<double>().compile(expression_string, expression);
}

Eigen::VectorXd Expression::evaluate(const Eigen::MatrixXd& r) {
    // define the output matrix
    Eigen::VectorXd output(r.rows());

    // loop over the vector of variables
    for (int i = 0; i < r.rows(); i++) {

        // set the variables to the expression
        for (int j = 0; j < r.cols(); j++) variables.at(j) = r(i, j);

        // evaluate the expression
        output(i) = expression.value();
    }

    // return the output matrix
    return output;
}

torch::Tensor Expression::evaluate(const torch::Tensor& r) {
    // define the output tensor and create the accessor
    torch::Tensor output = torch::zeros({r.size(0)}, torch::kDouble); torch::TensorAccessor<double, 1> output_accessor = output.accessor<double, 1>();

    // loop over the tensor of variables
    for (int i = 0; i < r.size(0); i++) {

        // set the variales to the expression
        std::memcpy(variables.data(), r.data_ptr<double>() + i * r.stride(0), r.stride(0) * sizeof(double));

        // evaluate the expression
        output_accessor[i] = expression.value();
    }

    // return the output tensor
    return output;
}
