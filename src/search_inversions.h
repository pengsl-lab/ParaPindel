#ifndef SEARCH_INVERSIONS_H
#define	SEARCH_INVERSIONS_H

int searchInversions(ControlState& currentState, unsigned NumBoxes, const SearchWindow& currentWindow, int mpirank);
void LeftMostINV(ControlState& currentState, SPLIT_READ & currentRead, const SearchWindow& window);
#endif /* SEARCH_INVERSIONS_H */
