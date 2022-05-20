# Comments from CRAN

still uninitialized use plus memory leaks in the package's own code
(and more in system libs which may reflect how they are used)

Please fix and resubmit.

# Response

Found the uninitialized value and retested everything using valgrind.

There were no uninitialized values on my system.  The following was the printout at the end:

==20283== LEAK SUMMARY:
==20283==    definitely lost: 0 bytes in 0 blocks
==20283==    indirectly lost: 0 bytes in 0 blocks
==20283==      possibly lost: 197,221 bytes in 6,400 blocks
==20283==    still reachable: 879,637,475 bytes in 212,314 blocks
==20283==                       of which reachable via heuristic:
==20283==                         newarray           : 4,264 bytes in 1 blocks
==20283==         suppressed: 0 bytes in 0 blocks
==20283== Reachable blocks (those to which a pointer was found) are not shown.
==20283== To see them, rerun with: --leak-check=full --show-leak-kinds=all
==20283== 
==20283== For lists of detected and suppressed errors, rerun with: -s
==20283== ERROR SUMMARY: 304 errors from 304 contexts (suppressed: 0 from 0)


All of these memory issues came from udunits but not anything from
within nlmixr2est.
