#include "Main.h"
#include "Declaration of procedures.h"

int main()
{

    Iter_Glob = 0;
    Create_WriteDir();
    Mesh_Init();
    Initial_Conditions();

    Section_value_MUSCL(xx_1, yy_1, "NULL");

    Redistricting1();
    Write();

    double test = 0.0;

    do
    {
        Redistricting();

        //for (int Iter = 0; Iter <= 100; Iter++)
        {

            //Test(Iter);

            Iter_Glob++;

            Redistricting1();
            Calculation_Velocity_U();
            Redistricting1();
            Calculation_Pressure_P();
            Redistricting1();

        }

        Development();
        Write();

        Time();
        
    } while (_time <= 0.0);

    Write_End();
}