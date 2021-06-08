/*
 * This File is part of Pindel; a program to locate genomic variation.
 * https://trac.nbic.nl/pindel/
 *
 *   Copyright (C) 2011 Kai Ye
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef SEARCHER_H
#define	SEARCHER_H

bool CheckMismatches (const std::string & TheInput, const std::string & CurrentReadSeq, const UniquePoint & UP, short & numberOfMismatch);
void CategorizePositions(const char readBase, const std::string & chromosomeSeq, const std::vector<PosVector>& PD_Plus, std::vector<PosVector>& PD_Plus_Output, const int numMisMatches,
                         const int searchDirection,	const int maxNumMismatches );

/*void CheckLeft_Far (SPLIT_READ & OneRead,
										const std::string & TheInput,
										const std::string & CurrentReadSeq,
										const std::vector < unsigned int >Left_PD[],
										const short &BP_Left_Start,
										const short &BP_Left_End,
										const short &CurrentLength,
										SortedUniquePoints &LeftUP);*/

/*void CheckRight_Far (SPLIT_READ & OneRead,
										 const std::string & TheInput,
										 const std::string & CurrentReadSeq,
										 const std::vector < unsigned int >Right_PD[],
										 const short &BP_Right_Start,
										 const short &BP_Right_End,
										 const short &CurrentPos,
										 SortedUniquePoints &RightUP);*/

void CheckLeft_Close (SPLIT_READ & OneRead,
                      const std::string & TheInput,
                      const std::string & CurrentReadSeq,
                      const std::vector < PosVector >& Left_PD,
                      const short &BP_Left_Start,
                      const short &BP_Left_End,
                      const short &CurrentLength,
                      SortedUniquePoints &LeftUP);

void CheckRight_Close (SPLIT_READ & OneRead,
                       const std::string & TheInput,
                       const std::string & CurrentReadSeq,
                       const std::vector < PosVector >& Right_PD,
                       const short &BP_Right_Start,
                       const short &BP_Right_End,
                       const short &CurrentPos,
                       SortedUniquePoints &RightUP);
void CheckLeft_Close_Perfect (SPLIT_READ & OneRead,
                              const std::string & TheInput,
                              const std::string & CurrentReadSeq,
                              const std::vector < PosVector >& Left_PD,
                              const short &BP_Left_Start,
                              const short &BP_Left_End,
                              const short &CurrentLength,
                              SortedUniquePoints &LeftUP);

void CheckRight_Close_Perfect (SPLIT_READ & OneRead,
                               const std::string & TheInput,
                               const std::string & CurrentReadSeq,
                               const std::vector < PosVector >& Right_PD,
                               const short &BP_Right_Start,
                               const short &BP_Right_End,
                               const short &CurrentPos,
                               SortedUniquePoints &RightUP);
#endif /* SEARCHER_H */
