#include <exprtk.hpp>

struct Expression {
    exprtk::expression<double>   expression;
    exprtk::symbol_table<double>    symbols;
    std::vector<double>                vars;
};

extern "C" {
    typedef unsigned long ulong;

    void* compile(const char* string, ulong nvars) {
        Expression* expr = new Expression(); exprtk::parser<double> parser; expr->vars.resize(nvars);

        for (int i = 1; i <= nvars; i++) {
            expr->symbols.add_variable("r" + std::to_string(i), expr->vars.at(i - 1));
        }

        expr->symbols.add_constants(); expr->expression.register_symbol_table(expr->symbols);

        if (!parser.compile(string, expr->expression)) {
            throw std::runtime_error("FAILED TO COMPILE EXPRESSION: " + std::string(parser.error().c_str()));
        }

        return static_cast<void*>(expr);
    }

    double evaluate(void* exprv, double* vars) {
        Expression* expr = static_cast<Expression*>(exprv);

        for (size_t i = 1; i <= expr->vars.size(); i++) {
            expr->symbols.get_variable("r" + std::to_string(i))->ref() = vars[i - 1];
        }

        return expr->expression.value();
    }

    void deinit(void* exprv) {
        delete static_cast<Expression*>(exprv);
    }
}
