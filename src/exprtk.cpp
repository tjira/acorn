#include <exprtk.hpp>

struct Expression {
    exprtk::expression<double>   expression;
    exprtk::symbol_table<double>    symbols;
    std::vector<double>                vars;
};

extern "C" {
    void* compile(const char* string, int nvars) {
        Expression* expr = new Expression(); exprtk::parser<double> parser; expr->vars.resize(nvars);

        for (int i = 1; i <= nvars; i++) {
            expr->symbols.add_variable("r" + std::to_string(i), expr->vars.at(i - 1));
        }

        expr->symbols.add_constants(); expr->expression.register_symbol_table(expr->symbols);

        parser.compile(string, expr->expression);

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
