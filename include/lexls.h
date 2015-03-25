/*
  API for internal::LexLSE and internal::LexLSI
*/

#ifndef LEXLSE_API
#define LEXLSE_API

namespace LexLS
{

    class LexLSE
    {
    public:

        LexLSE(){}

        LexLSE(Index nVar_, Index nObj_, Index *ObjDim_)
        {
            resize(nVar_, nObj_, ObjDim_);
            setObjDim(ObjDim_);
        }

        void resize(Index nVar_, Index nObj_, Index *ObjDim_)
        {
            lexlse.resize(nVar_, nObj_, ObjDim_);
        }

        void setObjDim(Index *ObjDim_)
        {
            lexlse.setObjDim(ObjDim_);
        }

        // @todo introduce enum
        // @todo include solve_option in the list of parameters
        dVectorType& solve(Index solve_option = 0)
        {
            lexlse.factorize();

            switch(solve_option){

            case 0:
                lexlse.solve();
                break;

            case 1:
                lexlse.solveLeastNorm_1();
                break;

            case 2:
                lexlse.solveLeastNorm_2();
                break;

            case 3:
                lexlse.solveLeastNorm_3();
                break;
            }

            return lexlse.get_x();
        }

    private:

        internal::LexLSE lexlse;
    };

}

#endif // LEXLSE_API
