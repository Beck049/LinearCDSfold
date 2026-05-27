#ifndef VIE_INF
#define VIE_INF 999999999
#endif
#ifndef NBPAIRS
#define NBPAIRS 7
#endif

/* Internal Loop 2x1 / 2x1 的 Free Energy
 * 1st NBPAIRS : ij 的 pair
 * 2nd NBPAIRS : pq 的 pair
 * 3nd, 4th, 5th : 對應 nuc
 * 
 * Usage:
 *  If 左側為 1
 *    int21_37[type][type_2][nuci1][nucq1][nucj_1];
 *  If 右側為 1
 *    int21_37[type_2][type][nucq1][nuci1][nucp_1];
 */ 

int int21_37[NBPAIRS+1][NBPAIRS+1][5][5][5] =
{{{{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  ,{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  ,{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  ,{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  ,{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  }
 ,{{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  ,{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  ,{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  ,{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  ,{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  }
 ,{{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  ,{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  ,{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  ,{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  ,{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  }
 ,{{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  ,{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  ,{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  ,{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  ,{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  }
 ,{{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  ,{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  ,{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  ,{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  ,{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  }
 ,{{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  ,{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  ,{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  ,{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  ,{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  }
 ,{{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  ,{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  ,{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  ,{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  ,{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  }
 ,{{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  ,{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  ,{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  ,{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  ,{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  }
 }
,{{{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  ,{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  ,{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  ,{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  ,{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  }
 ,{{{   230,   230,   230,   230,   230}
   ,{   230,   230,   230,   230,   230}
   ,{   230,   230,   230,   230,   230}
   ,{   230,   230,   230,   230,   230}
   ,{   230,   230,   230,   230,   230}
   }
  ,{{   230,   230,   230,   110,   230}
   ,{   230,   230,   230,   110,   230}
   ,{   230,   230,   230,   110,   230}
   ,{   110,   110,   110,   110,   110}
   ,{   230,   230,   230,   110,   230}
   }
  ,{{   230,   230,   230,   230,   230}
   ,{   230,   230,   230,   230,   230}
   ,{   230,   230,   230,   230,   230}
   ,{   230,   230,   230,   230,   230}
   ,{   230,   230,   230,   230,   230}
   }
  ,{{   230,   110,   230,   110,   230}
   ,{   110,   110,   110,   110,   110}
   ,{   230,   110,   230,   110,   230}
   ,{   110,   110,   110,   110,   110}
   ,{   230,   110,   230,   110,   230}
   }
  ,{{   230,   230,   230,   230,   150}
   ,{   230,   230,   230,   230,   150}
   ,{   230,   230,   230,   230,   150}
   ,{   230,   230,   230,   230,   150}
   ,{   150,   150,   150,   150,   150}
   }
  }
 ,{{{   250,   250,   250,   230,   230}
   ,{   250,   250,   230,   230,   230}
   ,{   250,   230,   250,   230,   230}
   ,{   230,   230,   230,   230,   230}
   ,{   250,   250,   230,   230,   230}
   }
  ,{{   250,   250,   230,   110,   230}
   ,{   250,   250,   230,   110,   230}
   ,{   230,   230,   170,   110,   230}
   ,{   110,    80,   110,   110,   110}
   ,{   230,   230,   230,   110,   230}
   }
  ,{{   250,   250,   250,   230,   230}
   ,{   230,   230,   230,   230,   230}
   ,{   250,   230,   250,   230,   230}
   ,{   230,   230,   230,   230,   230}
   ,{   250,   250,   230,   230,   230}
   }
  ,{{   230,   170,   230,   110,   230}
   ,{   230,   170,   230,    80,   230}
   ,{   230,   110,   230,   110,   230}
   ,{   120,   120,   110,   110,   110}
   ,{   230,   110,   230,   110,   230}
   }
  ,{{   230,   230,   230,   230,   150}
   ,{   230,   230,   230,   230,   150}
   ,{   230,   230,   220,   230,   150}
   ,{   230,   230,   230,   230,   150}
   ,{   170,   150,   170,   150,   140}
   }
  }
 ,{{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   300,   300,   190,   300}
   ,{   300,   300,   300,   190,   300}
   ,{   300,   300,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   300,   300,   190,   300}
   }
  ,{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   190,   300,   190,   300}
   ,{   300,   190,   300,   190,   300}
   ,{   300,   190,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   190,   300,   190,   300}
   }
  ,{{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   220,   220,   220,   220,   220}
   }
  }
 ,{{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   300,   300,   190,   300}
   ,{   300,   300,   300,   190,   300}
   ,{   300,   300,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   300,   300,   190,   300}
   }
  ,{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   190,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   190,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   190,   300,   190,   300}
   }
  ,{{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   220,   220,   220,   220,   220}
   }
  }
 ,{{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   300,   300,   190,   300}
   ,{   300,   300,   300,   190,   300}
   ,{   300,   300,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   300,   300,   190,   300}
   }
  ,{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   190,   300,   190,   300}
   ,{   300,   190,   300,   190,   300}
   ,{   300,   190,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   190,   300,   190,   300}
   }
  ,{{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   220,   220,   220,   220,   220}
   }
  }
 ,{{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   300,   300,   190,   300}
   ,{   300,   300,   300,   190,   300}
   ,{   300,   300,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   300,   300,   190,   300}
   }
  ,{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   190,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   190,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   190,   300,   190,   300}
   }
  ,{{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   220,   220,   220,   220,   220}
   }
  }
 ,{{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   300,   300,   190,   300}
   ,{   300,   300,   300,   190,   300}
   ,{   300,   300,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   300,   300,   190,   300}
   }
  ,{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   190,   300,   190,   300}
   ,{   300,   190,   300,   190,   300}
   ,{   300,   190,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   190,   300,   190,   300}
   }
  ,{{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   220,   220,   220,   220,   220}
   }
  }
 }
,{{{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  ,{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  ,{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  ,{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  ,{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  }
 ,{{{   250,   250,   230,   230,   230}
   ,{   250,   250,   230,   230,   230}
   ,{   230,   230,   230,   230,   230}
   ,{   230,   230,   230,   230,   230}
   ,{   230,   230,   230,   230,   230}
   }
  ,{{   250,   250,   230,   230,   230}
   ,{   250,   250,   230,   210,   230}
   ,{   230,   230,   230,   230,   230}
   ,{   120,   120,   110,   110,   110}
   ,{   230,   230,   230,   230,   230}
   }
  ,{{   230,   230,   230,   230,   230}
   ,{   230,   230,   230,   230,   230}
   ,{   230,   230,   230,   230,   230}
   ,{   230,   230,   230,   230,   230}
   ,{   230,   230,   190,   230,   230}
   }
  ,{{   230,   110,   230,   110,   230}
   ,{   110,   110,   110,   110,   110}
   ,{   230,   110,   230,   110,   230}
   ,{   110,   110,   110,   110,   110}
   ,{   230,   110,   230,   110,   230}
   }
  ,{{   230,   230,   230,   230,   150}
   ,{   230,   230,   230,   230,   150}
   ,{   230,   230,   230,   230,   150}
   ,{   230,   230,   230,   230,   150}
   ,{   150,   150,   150,   150,   150}
   }
  }
 ,{{{   230,   230,   230,   230,   230}
   ,{   230,   230,   230,   230,   230}
   ,{   230,   230,   230,   230,   230}
   ,{   230,   230,   230,   230,   230}
   ,{   230,   230,   230,   230,   230}
   }
  ,{{   230,   230,   230,   230,   230}
   ,{   230,   230,   230,   230,   230}
   ,{   230,   230,   230,   230,   230}
   ,{   110,   110,   110,   110,   110}
   ,{   230,   230,   230,   230,   230}
   }
  ,{{   230,   230,   230,   230,   230}
   ,{   230,   230,   230,   230,   230}
   ,{   230,   230,   230,   230,   230}
   ,{   230,   230,   230,   230,   230}
   ,{   230,   230,   230,   230,   230}
   }
  ,{{   230,   110,   230,   110,   230}
   ,{   230,   110,   230,   110,   230}
   ,{   230,   110,   230,   110,   230}
   ,{   110,   110,   110,   110,   110}
   ,{   230,   110,   230,   110,   230}
   }
  ,{{   230,   230,   230,   230,   150}
   ,{   230,   230,   230,   230,   150}
   ,{   230,   230,   230,   230,   150}
   ,{   230,   230,   230,   230,   150}
   ,{   150,   150,   150,   150,   150}
   }
  }
 ,{{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   190,   300,   190,   300}
   ,{   300,   190,   300,   190,   300}
   ,{   300,   190,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   190,   300,   190,   300}
   }
  ,{{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   220,   220,   220,   220,   220}
   }
  }
 ,{{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   300,   300,   300,   300}
   ,{   300,   250,   300,   210,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   190,   120,   190,   190,   190}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   190,   300,   300}
   }
  ,{{   300,   190,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   190,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   190,   300,   190,   300}
   }
  ,{{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   220,   220,   220,   220,   220}
   }
  }
 ,{{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   190,   300,   190,   300}
   ,{   300,   190,   300,   190,   300}
   ,{   300,   190,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   190,   300,   190,   300}
   }
  ,{{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   220,   220,   220,   220,   220}
   }
  }
 ,{{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   190,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   190,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   190,   300,   190,   300}
   }
  ,{{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   220,   220,   220,   220,   220}
   }
  }
 ,{{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   190,   300,   190,   300}
   ,{   300,   190,   300,   190,   300}
   ,{   300,   190,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   190,   300,   190,   300}
   }
  ,{{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   220,   220,   220,   220,   220}
   }
  }
 }
,{{{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  ,{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  ,{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  ,{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  ,{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  }
 ,{{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   300,   300,   300,   300}
   ,{   300,   250,   300,   210,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   190,   120,   190,   190,   190}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   190,   300,   300}
   }
  ,{{   300,   190,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   190,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   190,   300,   190,   300}
   }
  ,{{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   220,   220,   220,   220,   220}
   }
  }
 ,{{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   190,   300,   190,   300}
   ,{   300,   190,   300,   190,   300}
   ,{   300,   190,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   190,   300,   190,   300}
   }
  ,{{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   220,   220,   220,   220,   220}
   }
  }
 ,{{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   260,   370,   260,   370}
   ,{   370,   260,   370,   260,   370}
   ,{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   300,   300,   300,   300,   300}
   }
  }
 ,{{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   250,   370,   210,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   260,   120,   260,   260,   260}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   190,   370,   370}
   }
  ,{{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   300,   300,   300,   300,   300}
   }
  }
 ,{{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   260,   370,   260,   370}
   ,{   370,   260,   370,   260,   370}
   ,{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   300,   300,   300,   300,   300}
   }
  }
 ,{{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   300,   300,   300,   300,   300}
   }
  }
 ,{{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   260,   370,   260,   370}
   ,{   370,   260,   370,   260,   370}
   ,{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   300,   300,   300,   300,   300}
   }
  }
 }
,{{{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  ,{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  ,{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  ,{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  ,{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  }
 ,{{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   300,   300,   190,   300}
   ,{   300,   300,   300,   190,   300}
   ,{   300,   300,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   300,   300,   190,   300}
   }
  ,{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   190,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   190,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   190,   300,   190,   300}
   }
  ,{{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   220,   220,   220,   220,   220}
   }
  }
 ,{{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   300,   300,   190,   300}
   ,{   300,   300,   300,   190,   300}
   ,{   300,   300,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   300,   300,   190,   300}
   }
  ,{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   190,   300,   190,   300}
   ,{   300,   190,   300,   190,   300}
   ,{   300,   190,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   190,   300,   190,   300}
   }
  ,{{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   220,   220,   220,   220,   220}
   }
  }
 ,{{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   260,   370}
   ,{   370,   370,   370,   260,   370}
   ,{   370,   370,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   370,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   260,   370,   260,   370}
   ,{   370,   260,   370,   260,   370}
   ,{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   300,   300,   300,   300,   300}
   }
  }
 ,{{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   260,   370}
   ,{   370,   370,   370,   260,   370}
   ,{   370,   370,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   370,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   300,   300,   300,   300,   300}
   }
  }
 ,{{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   260,   370}
   ,{   370,   370,   370,   260,   370}
   ,{   370,   370,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   370,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   260,   370,   260,   370}
   ,{   370,   260,   370,   260,   370}
   ,{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   300,   300,   300,   300,   300}
   }
  }
 ,{{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   260,   370}
   ,{   370,   370,   370,   260,   370}
   ,{   370,   370,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   370,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   300,   300,   300,   300,   300}
   }
  }
 ,{{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   260,   370}
   ,{   370,   370,   370,   260,   370}
   ,{   370,   370,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   370,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   260,   370,   260,   370}
   ,{   370,   260,   370,   260,   370}
   ,{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   300,   300,   300,   300,   300}
   }
  }
 }
,{{{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  ,{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  ,{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  ,{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  ,{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  }
 ,{{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   190,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   190,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   190,   300,   190,   300}
   }
  ,{{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   220,   220,   220,   220,   220}
   }
  }
 ,{{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   190,   300,   190,   300}
   ,{   300,   190,   300,   190,   300}
   ,{   300,   190,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   190,   300,   190,   300}
   }
  ,{{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   220,   220,   220,   220,   220}
   }
  }
 ,{{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   260,   370,   260,   370}
   ,{   370,   260,   370,   260,   370}
   ,{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   300,   300,   300,   300,   300}
   }
  }
 ,{{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   300,   300,   300,   300,   300}
   }
  }
 ,{{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   260,   370,   260,   370}
   ,{   370,   260,   370,   260,   370}
   ,{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   300,   300,   300,   300,   300}
   }
  }
 ,{{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   300,   300,   300,   300,   300}
   }
  }
 ,{{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   260,   370,   260,   370}
   ,{   370,   260,   370,   260,   370}
   ,{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   300,   300,   300,   300,   300}
   }
  }
 }
,{{{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  ,{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  ,{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  ,{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  ,{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  }
 ,{{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   300,   300,   190,   300}
   ,{   300,   300,   300,   190,   300}
   ,{   300,   300,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   300,   300,   190,   300}
   }
  ,{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   190,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   190,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   190,   300,   190,   300}
   }
  ,{{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   220,   220,   220,   220,   220}
   }
  }
 ,{{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   300,   300,   190,   300}
   ,{   300,   300,   300,   190,   300}
   ,{   300,   300,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   300,   300,   190,   300}
   }
  ,{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   190,   300,   190,   300}
   ,{   300,   190,   300,   190,   300}
   ,{   300,   190,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   190,   300,   190,   300}
   }
  ,{{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   220,   220,   220,   220,   220}
   }
  }
 ,{{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   260,   370}
   ,{   370,   370,   370,   260,   370}
   ,{   370,   370,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   370,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   260,   370,   260,   370}
   ,{   370,   260,   370,   260,   370}
   ,{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   300,   300,   300,   300,   300}
   }
  }
 ,{{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   260,   370}
   ,{   370,   370,   370,   260,   370}
   ,{   370,   370,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   370,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   300,   300,   300,   300,   300}
   }
  }
 ,{{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   260,   370}
   ,{   370,   370,   370,   260,   370}
   ,{   370,   370,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   370,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   260,   370,   260,   370}
   ,{   370,   260,   370,   260,   370}
   ,{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   300,   300,   300,   300,   300}
   }
  }
 ,{{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   260,   370}
   ,{   370,   370,   370,   260,   370}
   ,{   370,   370,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   370,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   300,   300,   300,   300,   300}
   }
  }
 ,{{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   260,   370}
   ,{   370,   370,   370,   260,   370}
   ,{   370,   370,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   370,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   260,   370,   260,   370}
   ,{   370,   260,   370,   260,   370}
   ,{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   300,   300,   300,   300,   300}
   }
  }
 }
,{{{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  ,{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  ,{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  ,{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  ,{{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   ,{   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF,   VIE_INF}
   }
  }
 ,{{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   190,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   190,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   190,   300,   190,   300}
   }
  ,{{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   220,   220,   220,   220,   220}
   }
  }
 ,{{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   190,   300,   190,   300}
   ,{   300,   190,   300,   190,   300}
   ,{   300,   190,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   190,   300,   190,   300}
   }
  ,{{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   220,   220,   220,   220,   220}
   }
  }
 ,{{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   260,   370,   260,   370}
   ,{   370,   260,   370,   260,   370}
   ,{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   300,   300,   300,   300,   300}
   }
  }
 ,{{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   300,   300,   300,   300,   300}
   }
  }
 ,{{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   260,   370,   260,   370}
   ,{   370,   260,   370,   260,   370}
   ,{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   300,   300,   300,   300,   300}
   }
  }
 ,{{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   300,   300,   300,   300,   300}
   }
  }
 ,{{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   260,   370,   260,   370}
   ,{   370,   260,   370,   260,   370}
   ,{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   300,   300,   300,   300,   300}
   }
  }
 }};
