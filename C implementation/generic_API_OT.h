
void print_bytes(unsigned char *p, size_t s);

int generate_2n_basis(const f2elm_t, const f2elm_t, const f2elm_t,
                      point_proj_t, point_proj_t,
                      const f2elm_t, const f2elm_t, const f2elm_t);
            
/*                      
int generate_2n_basis(const f2elm_t*, const f2elm_t*, const f2elm_t*,
                      point_proj_t*, point_proj_t*,
                      const f2elm_t*, const f2elm_t*, const f2elm_t*);
*/
                      
void pointAddition(const point_proj*, const point_proj*,
                   const f2elm_t*, const f2elm_t*, const f2elm_t*,
                   point_proj*);
                   
void pointDifference(const point_proj*, const point_proj*,
                     const f2elm_t*, const f2elm_t*, const f2elm_t*,
                     point_proj*);
                
void projectiveSum(const point_proj*, const point_proj*, const f2elm_t*, point_full_proj*);

void getProjective_Y(const point_proj*, const f2elm_t*, f2elm_t*);

