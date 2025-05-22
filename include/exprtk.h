typedef unsigned long ulong;

void*   compile(const char* string, ulong   nvars);
double evaluate(void*        exprv, double* vars );
void     deinit(void*        exprv               );
