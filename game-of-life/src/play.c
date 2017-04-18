#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <omp.h>
#include <game-of-life.h>
#include <mpi.h>

void play (int *board, int *newboard, int N,int M) {
	/*
	(copied this from some web page, hence the English spellings...)

	1.STASIS : If, for a given cell, the number of on neighbours is
	exactly two, the cell maintains its status quo into the next
	generation. If the cell is on, it stays on, if it is off, it stays off.

	2.GROWTH : If the number of on neighbours is exactly three, the cell
	will be on in the next generation. This is regardless of the cell's
	current state.

	3.DEATH : If the number of on neighbours is 0, 1, 4-8, the cell will
	be off in the next generation.
	*/
	struct timeval startwtime, endwtime, startwtimetemp, endwtimetemp; // time counting variables
	double rearrange_time;
	MPI_Datatype rowtype;
	MPI_Type_contiguous(M, MPI_INT, &rowtype);
	MPI_Type_commit(&rowtype);

	gettimeofday (&startwtime, NULL);//start timer

	int   i, j, a,numtasks,rank;
	MPI_Comm_size(MPI_COMM_WORLD,&numtasks); //get number of tasks
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);//get my rank

	/* for each cell, apply the rules of Life */
	if (numtasks==1)
	{
		#pragma omp parallel for private(i,j,a)
		for (i=0; i<N; i++)
		for (j=0; j<M; j++) {
			a = adjacent_to (board, i, j, N,M);
			if (a == 2) NewBoard(i,j) = Board(i,j);
			if (a == 3) NewBoard(i,j) = 1;
			if (a < 2) NewBoard(i,j) = 0;
			if (a > 3) NewBoard(i,j) = 0;
		}

	}
	else if (numtasks==2)
	{
		MPI_Request reqs[4];   // required variable for non-blocking calls
		MPI_Status stats[4];
		int UP[M], DOWN[M];


		MPI_Isend(&board[(N-1)*M], 1, rowtype, (rank+1)%2 , 2, MPI_COMM_WORLD, &reqs[1]);
		MPI_Isend(&board[0], 1, rowtype, (rank+1)%2 , 1, MPI_COMM_WORLD, &reqs[0]);

		MPI_Irecv(&UP[0], M, MPI_INT, (rank+1)%2, 2, MPI_COMM_WORLD, &reqs[2]);
		MPI_Irecv(&DOWN[0], M, MPI_INT, (rank+1)%2, 1, MPI_COMM_WORLD, &reqs[3]);


		gettimeofday (&startwtimetemp, NULL);
		#pragma omp parallel for private(i,j,a)
		for (i=1; i<N-1; i++)  //h prwth kai h teleytaia grammh den xreiazetai mpi giati yparxoyn ola ta dedomena topika
		for (j=0; j<M; j++) {
			a = adjacent_to (board, i, j, N,M);
			if (a == 2) NewBoard(i,j) = Board(i,j);
			if (a == 3) NewBoard(i,j) = 1;
			if (a < 2) NewBoard(i,j) = 0;
			if (a > 3) NewBoard(i,j) = 0;
		}

		MPI_Waitall(4, reqs, stats);
		//twra mporw na parw thn prwth kai thn teleytaia grammh

		//kanontas thn 1h grammh
		#pragma omp parallel for private(j,a)
		for (j=0 ; j<M; j++)
		{
			a = 0;
			if (j==0)
			{
				a = a + UP[0];
				a = a + UP[M-1];
				a = a + UP[1];
				a = a + Board(0,M-1);
				a = a + Board(0,1);
				a = a + Board(1,M-1);
				a = a + Board(1,0);
				a = a + Board(1,1);
			}else if(j == (M-1))
			{
				a = a + UP[M-1];
				a = a + UP[M-2];
				a = a + UP[0];
				a = a + Board(0,M-2);
				a = a + Board(0,0);
				a = a + Board(1,M-1);
				a = a + Board(1,M-2);
				a = a + Board(1,0);

			}else
			{
				//  board(0,j) //to stoixeio maintains
				a = a + Board(0,j-1);
				a = a + Board(0,j+1);
				a = a + Board(1,j-1);
				a = a + Board(1,j);
				a = a + Board(1,j+1);
				a = a + UP[j-1];
				a = a + UP[j];
				a = a + UP[j+1];
			}
			//analoga me toys geitones moy zw h pethainw
			if (a == 2) NewBoard(0,j) = Board(0,j);
			if (a == 3) NewBoard(0,j) = 1;
			if (a < 2) NewBoard(0,j) = 0;
			if (a > 3) NewBoard(0,j) = 0;
		}

		//kanontas thn grammh N-1
		#pragma omp parallel for private(j,a)
		for (j=0 ; j<M; j++)
		{
			a = 0;
			if (j==0)
			{
				a = a + DOWN[0];
				a = a + DOWN[M-1];
				a = a + DOWN[1];
				a = a + Board(N-1,M-1);
				a = a + Board(N-1,1);
				a = a + Board(N-2,M-1);
				a = a + Board(N-2,0);
				a = a + Board(N-2,1);
			}else if(j == M-1)
			{
				a = a + DOWN[M-1];
				a = a + DOWN[M-2];
				a = a + DOWN[0];
				a = a + Board(N-1,M-2);
				a = a + Board(N-1,0);
				a = a + Board(N-2,M-1);
				a = a + Board(N-2,M-2);
				a = a + Board(N-2,0);

			}else
			{
				a = a + Board(N-1,j-1);
				a = a + Board(N-1,j+1);
				a = a + Board(N-2,j-1);
				a = a + Board(N-2,j);
				a = a + Board(N-2,j+1);
				a = a + DOWN[j-1];
				a = a + DOWN[j];
				a = a + DOWN[j+1];
			}
			//analoga me toys geitones moy zw h pethainw
			if (a == 2) NewBoard(N-1,j) = Board(N-1,j);
			if (a == 3) NewBoard(N-1,j) = 1;
			if (a < 2) NewBoard(N-1,j) = 0;
			if (a > 3) NewBoard(N-1,j) = 0;
		}


	}else if(numtasks==4)
	{
		gettimeofday (&startwtime, NULL);// tld

		int sourceU, sourceD, destU, destD;  //milame gia to sygkekrimmeno task
		MPI_Request reqs[4];   // required variable for non-blocking calls
		MPI_Status stats[4];
		int UP[M], DOWN[M];

		if (rank==0)
		{
			sourceU = 3;
			sourceD = 1;
			destU = 3;
			destD = 1;
		}else if(rank ==1)
		{
			sourceU = 0;
			sourceD = 2;
			destU = 0;
			destD = 2;

		}else if(rank==2)
		{
			sourceU = 1;
			sourceD = 3;
			destU = 1;
			destD = 3;
		}else if(rank==3)
		{
			sourceU = 2;
			sourceD = 0;
			destU = 2;
			destD = 0;
		}


		MPI_Isend(&board[(N-1)*M], 1, rowtype, destD, 1, MPI_COMM_WORLD, &reqs[1]);  //any tag giati dexetai appo diafotika processes
		MPI_Isend(&board[0], 1, rowtype, destU, 1, MPI_COMM_WORLD, &reqs[0]);

		MPI_Irecv(&UP[0], M, MPI_INT, sourceU, 1, MPI_COMM_WORLD, &reqs[2]);
		MPI_Irecv(&DOWN[0], M, MPI_INT, sourceD, 1, MPI_COMM_WORLD, &reqs[3]);

			gettimeofday (&startwtimetemp, NULL);
		#pragma omp parallel for private(i,j,a)
		for (i=1; i<N-1; i++)  //h prwth kai h teleytaia grammh den xreiazetai mpi giati yparxoyn ola ta dedomena topika
		for (j=0; j<M; j++) {
			a = adjacent_to (board, i, j, N,M);
			if (a == 2) NewBoard(i,j) = Board(i,j);
			if (a == 3) NewBoard(i,j) = 1;
			if (a < 2) NewBoard(i,j) = 0;
			if (a > 3) NewBoard(i,j) = 0;
		}
		//twra mporw na parw thn prwth kai thn teleytaia grammh
		MPI_Waitall(4, reqs, stats);
		#pragma omp parallel for private(j,a)
		for (j=0 ; j<M; j++)
		{
			a = 0;
			if (j==0)
			{
				a = a + UP[0];
				a = a + UP[M-1];
				a = a + UP[1];
				a = a + Board(0,M-1);
				a = a + Board(0,1);
				a = a + Board(1,M-1);
				a = a + Board(1,0);
				a = a + Board(1,1);
			}else if(j == (M-1))
			{
				a = a + UP[M-1];
				a = a + UP[M-2];
				a = a + UP[0];
				a = a + Board(0,M-2);
				a = a + Board(0,0);
				a = a + Board(1,M-1);
				a = a + Board(1,M-2);
				a = a + Board(1,0);

			}else
			{
				//  board(0,j) //to stoixeio maintains
				a = a + Board(0,j-1);
				a = a + Board(0,j+1);
				a = a + Board(1,j-1);
				a = a + Board(1,j);
				a = a + Board(1,j+1);
				a = a + UP[j-1];
				a = a + UP[j];
				a = a + UP[j+1];
			}
			//analoga me toys geitones moy zw h pethainw
			if (a == 2) NewBoard(0,j) = Board(0,j);
			if (a == 3) NewBoard(0,j) = 1;
			if (a < 2) NewBoard(0,j) = 0;
			if (a > 3) NewBoard(0,j) = 0;
		}
		//kanontas thn grammh N-1
		#pragma omp parallel for private(j,a)
		for (j=0 ; j<M; j++)
		{
			a = 0;
			if (j==0)
			{
				a = a + DOWN[0];
				a = a + DOWN[M-1];
				a = a + DOWN[1];
				a = a + Board(N-1,M-1);
				a = a + Board(N-1,1);
				a = a + Board(N-2,M-1);
				a = a + Board(N-2,0);
				a = a + Board(N-2,1);
			}else if(j == M-1)
			{
				a = a + DOWN[M-1];
				a = a + DOWN[M-2];
				a = a + DOWN[0];
				a = a + Board(N-1,M-2);
				a = a + Board(N-1,0);
				a = a + Board(N-2,M-1);
				a = a + Board(N-2,M-2);
				a = a + Board(N-2,0);

			}else
			{
				a = a + Board(N-1,j-1);
				a = a + Board(N-1,j+1);
				a = a + Board(N-2,j-1);
				a = a + Board(N-2,j);
				a = a + Board(N-2,j+1);
				a = a + DOWN[j-1];
				a = a + DOWN[j];
				a = a + DOWN[j+1];
			}
			//analoga me toys geitones moy zw h pethainw
			if (a == 2) NewBoard(N-1,j) = Board(N-1,j);
			if (a == 3) NewBoard(N-1,j) = 1;
			if (a < 2) NewBoard(N-1,j) = 0;
			if (a > 3) NewBoard(N-1,j) = 0;

		}
	}
	// copy the new board back into the old board
	#pragma omp parallel for private(i,j)
	for (i=0; i<N; i++)
	for (j=0; j<M; j++) {
		Board(i,j) = NewBoard(i,j);
	}
	gettimeofday (&endwtime, NULL); //stop timer
	rearrange_time = (double)((endwtime.tv_usec - startwtime.tv_usec)
	/1.0e6 + endwtime.tv_sec - startwtime.tv_sec);
	printf("Time to play at task %d : %fs\n",rank,rearrange_time);


}
