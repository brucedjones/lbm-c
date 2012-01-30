#ifndef CGNS_OUTPUT_HANDLER_H
#define CGNS_OUTPUT_HANDLER_H

void cgns_open(void);
void cgns_write_grid(void);
void cgns_append_sol(int iter);
void cgns_append_sol_field(double field[9][17][21], char *name, int index_field);
void cgns_close(void);
void cgns_error_check(int error_code);

#endif