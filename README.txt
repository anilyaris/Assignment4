# ME-C3100 Computer Graphics, Fall 2017
# Lehtinen / Aarnio, Kemppinen, Ollikainen
#
# Assignment 4: Physical Simulation

Student name: Anıl YARIŞ
Student number: 643454
Hours spent on requirements (approx.): 4 hours
Hours spent on extra credit (approx.): 5 hours

# First, a 10-second poll about this assignment period:

Did you go to exercise sessions?
No.
Did you work on the assignment using Aalto computers, your own computers, or both?
Only on my computer.
# Which parts of the assignment did you complete? Mark them 'done'.
# You can also mark non-completed parts as 'attempted' if you spent a fair amount of
# effort on them. If you do, explain the work you did in the problems/bugs section
# and leave your 'attempt' code in place (commented out if necessary) so we can see it.

     R1 Euler integrator (1p): done
        R2 Spring system (2p): done
 R3 Trapezoid integrator (2p): done
      R4 Pendulum system (2p): done
         R5 Cloth system (3p): done

# Did you do any extra credit work?

RK4 implementation, particle generator (with bounce), random wind, cloth tearing.

Particle generator can be accessed from the sidebar similar to switching between other modes. The particles
generated are (unintentionally) in the form of "snakes" and they bounce off of the walls of a cube box whose 
dimensions can be adjusted with a variable within the code.

Cloth tearing can also be toggled from the sidebar. In my implementation the springs attached to the fixed
corners are not torn because otherwise they are always the first to tear and I thought this kind of defeats 
the purpose of cloth tearing.

For all of the randomly generated stuff the srand funcion does not appear within the particle_systems file so
it might seem like the rand() calls are not truly random, but in fact srand(time(NULL)) is called within the
constructor of App class. I implemented it this way so that the seed is set only once.

# Are there any known problems/bugs remaining in your code?

My simulations are not so stable. Even if I use the pre-implemented midpoint integration my systems explode
way faster compared to the example executable. This doesn't change much using RK4 either. My thought is I 
must have done something wrong in how the reaction forces should be reflected to other objects. I tried hard 
to figure out how the forces should act but I might not have found the correct way.

# Did you collaborate with anyone in the class?

No.

# Any other comments you'd like to share about the assignment or the course so far?

(Was the assignment too long? Too hard? Fun or boring? Did you learn something, or was it a total waste of time? Can we do something differently to help you learn? Please be brutally honest; we won't take it personally.)

