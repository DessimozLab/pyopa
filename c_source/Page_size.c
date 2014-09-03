/*
 * Page_size.c
 *
 *  Created on: Aug 20, 2014
 *      Author: machine
 */

#include "Page_size.h"

#if (defined _WIN32 || defined __WIN32__)
	#include <windows.h>
#else
	#include <unistd.h>
#endif

int getPageSize(void) {
#if (defined _WIN32 || defined __WIN32__)
	SYSTEM_INFO system_info;
	GetSystemInfo (&system_info);
	return system_info.dwPageSize;
#else
	return sysconf(_SC_PAGESIZE);
#endif
}
