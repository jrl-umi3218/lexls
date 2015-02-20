/*
  API for LexLSE
*/

#ifndef LEXLSE_API
#define LEXLSE_API

namespace LexLS
{

    class LexLSE
    {
    public:
        
        dVectorType& get_x()
        {
            return lexlse.get_x();
        }
        
    private:
        
        internal::LexLSE lexlse;
    };

}

#endif // LEXLSE_API
