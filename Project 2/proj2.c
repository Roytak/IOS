#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <signal.h>
#include <time.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <sys/wait.h>
#include <semaphore.h>

#define SIZE_OF_BLOCK 4096

int num_h, num_ox, max_queue_wait_time, max_molecule_wait_time;
FILE *shared_file;

/* shared variables */

int *n_action = NULL;
int *n_molecule = NULL;
int *n_h = NULL;
int *n_ox = NULL;
int *barrier_cnt = NULL;
int *h_cnt = NULL;
int *ox_cnt = NULL;

/* semaphores */

sem_t *write_file = NULL;
sem_t *h_q = NULL;
sem_t *ox_q = NULL;
sem_t *forming_h2o = NULL;
/* barrier */
sem_t *mutex = NULL;
sem_t *reu_bar_turnstile1 = NULL;
sem_t *reu_bar_turnstile2 = NULL;

int
parseArgs(int argc, char *argv[], int *num_h, int *num_ox, int *max_queue_wait_time, int *max_molecule_wait_time)
{
	const char *usage = "Usage : ./proj2 NO NH TI TB\n";
	char *ptr;
	long int arg1, arg2, arg3, arg4;

	if (argc != 5) {
		goto invalid_args;
	}

	for (int i = 1; i < argc; i++) {
		if (argv[i][0] == '\0') {
			goto invalid_args;
		}
		for (size_t j = 0; j < strlen(argv[i]); j++) {
			if (!isdigit(argv[i][j])) {
				goto invalid_args;
			}
		}
	}

	arg1 = strtoul(argv[1], &ptr, 10);
	arg2 = strtoul(argv[2], &ptr, 10);
	arg3 = strtoul(argv[3], &ptr, 10);
	arg4 = strtoul(argv[4], &ptr, 10);

	if (arg1 <= 0 || arg2 <= 0) {
		goto invalid_args;
	}

	if (arg3 < 0 || arg3 > 1000 || arg4 < 0 || arg4 > 1000) {
		goto invalid_args;
	}

	*num_ox = (int) arg1;
	*num_h = (int) arg2;
	*max_queue_wait_time = (int) arg3;
	*max_molecule_wait_time = (int) arg4;
	return 0;

	invalid_args:
		fprintf(stderr, "%s", usage);
		return 1;
}


/*
 *	creates a new mapping in the virtual address space and returns a pointer to it	
 */

void *
mapSharedMem(size_t size)
{
	void *ptr;

	/* let the kernel choose the place of the shared memory address by passing NULL as the first parameter */
	ptr = mmap(NULL, size, PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0); 

	if (ptr == (void *) -1) {
		perror(NULL);
		return NULL;
	}

	return ptr;
}

/*
 *	deletes the mappings on the address
 */

int
unmapSharedMem(void *mem_adress, size_t size)
{
	if (munmap(mem_adress, size)) {
		perror("Unmap");
		return -1;
	}
	return 0;
}

/*
 *	2 help functions for creating a barrier by using 3 semaphores
 */

void 
barrier_start()
{
	sem_wait(mutex);
	(*barrier_cnt) += 1;
	if (*barrier_cnt == 3) {
		sem_wait(reu_bar_turnstile2);
		sem_post(reu_bar_turnstile1);
	}
	sem_post(mutex);
	sem_wait(reu_bar_turnstile1);
	sem_post(reu_bar_turnstile1);
}

void
barrier_end()
{
	sem_wait(mutex);
	(*barrier_cnt) -= 1;
	if (*barrier_cnt == 0) {
		sem_wait(reu_bar_turnstile1);
		sem_post(reu_bar_turnstile2);
	}
	sem_post(mutex);
	sem_wait(reu_bar_turnstile2);
	sem_post(reu_bar_turnstile2);
}

/*
 * an oxygen process' function that writes into a shared file about its actions
 */

void
oxygen(int ox_id)
{
	int cur_ox_id = ox_id;
	unsigned int queue_wait_time = (rand() % (max_queue_wait_time + 1)) * 1000;
	unsigned int bonding_wait_time = (rand() % (max_molecule_wait_time + 1)) * 1000;

	sem_wait(write_file);
	fprintf(shared_file,"%d: O %d: started\n", (*n_action), ox_id);
	fflush(shared_file);
	(*n_action)++;
	sem_post(write_file);

	/* random backoff time after starting the process */
	usleep(queue_wait_time);

	sem_wait(write_file);
	fprintf(shared_file,"%d: O %d: going to queue\n", (*n_action), cur_ox_id);
	fflush(shared_file);
	(*n_action)++;
	sem_post(write_file);

	/* main mutex semaphore, it will not get released, until either a molecule is created or there is not enough hydrogens */
	sem_wait(forming_h2o);

	/* ran out of hydrogens to bond with */
	if ((*h_cnt + *n_h) <= 2) {
		sem_wait(write_file);
		fprintf(shared_file,"%d: O %d: not enough H\n", (*n_action), cur_ox_id);
		fflush(shared_file);
		(*n_action)++;
		sem_post(write_file);
		sem_post(forming_h2o);
		exit(0);		
	}

	(*n_ox) += 1;
	
	/* found hydrogens to bond with */
	if (*n_h >= 2) {
		sem_post(h_q);
		sem_post(h_q);
		(*n_h) -= 2;
		(*h_cnt) -= 2;
		sem_post(ox_q);
		(*n_ox) -= 1;
		(*ox_cnt)--;
	}
	else {
		sem_post(forming_h2o);
	}	


	sem_wait(ox_q);
	(*n_molecule)++;
	barrier_start();

	/* bonding */
	sem_wait(write_file);
	fprintf(shared_file,"%d: O %d: creating molecule %d\n", (*n_action), cur_ox_id, (*n_molecule));
	fflush(shared_file);
	(*n_action)++;
	sem_post(write_file);

	/* random backoff time for simulating the bonding */
	usleep(bonding_wait_time);

	barrier_end();

	sem_wait(write_file);
	fprintf(shared_file,"%d: O %d: molecule %d created\n", (*n_action), cur_ox_id, (*n_molecule));
	fflush(shared_file);
	(*n_action)++;
	sem_post(write_file);

	/* let oxygen release the main mutex */
	sem_post(forming_h2o);
	exit(0);
}

/*
 *	a hydrogen process' function
 */

void
hydrogen(int h_id)
{
	int cur_h_id = h_id;
	unsigned int queue_wait_time = (rand() % (max_queue_wait_time + 1)) * 1000;

	sem_wait(write_file);
	fprintf(shared_file,"%d: H %d: started\n", (*n_action),cur_h_id);
	fflush(shared_file);
	(*n_action)++;
	sem_post(write_file);

	/* random backoff time after starting the process */
	usleep(queue_wait_time);

	sem_wait(write_file);
	fprintf(shared_file,"%d: H %d: going to queue\n", (*n_action), cur_h_id);
	fflush(shared_file);
	(*n_action)++;
	sem_post(write_file);

	/* same mutex as in the oxygen process, the only difference being that hydrogens will not open it after bonding */
	sem_wait(forming_h2o);

	sem_wait(write_file);	
	sem_post(write_file);

	if ((*h_cnt) < 2 || (*ox_cnt) < 1) {
		sem_wait(write_file);
		fprintf(shared_file,"%d: H %d: not enough O or H\n", (*n_action), cur_h_id);
		fflush(shared_file);
		(*n_action)++;
		sem_post(write_file);
		sem_post(forming_h2o);
		exit(0);		
	} 

	(*n_h) += 1;

	/* found a hydrogen and an oxygen */
	if (*n_h >= 2 && *n_ox >= 1) {
		sem_post(h_q);
		sem_post(h_q);
		(*n_h) -= 2;
		(*h_cnt) -= 2;
		sem_post(ox_q);
		(*n_ox) -= 1;
	}
	else {
		sem_post(forming_h2o);
	}

	sem_wait(h_q);

	barrier_start();

	sem_wait(write_file);
	fprintf(shared_file,"%d: H %d: creating molecule %d\n", (*n_action), cur_h_id, (*n_molecule));
	fflush(shared_file);
	(*n_action)++;
	sem_post(write_file);

	barrier_end();

	sem_wait(write_file);
	fprintf(shared_file,"%d: H %d: molecule %d created\n", (*n_action), cur_h_id, (*n_molecule));
	fflush(shared_file);
	(*n_action)++;
	sem_post(write_file);

	exit(0);
}

/*
 *	unmap the sared memory and properly close all semaphores
 */

void
cleanup()
{
	unmapSharedMem((void *) n_action, sizeof(*n_action));
	unmapSharedMem((void *) n_molecule, sizeof(*n_molecule));
	unmapSharedMem((void *) n_h, sizeof(*n_h));
	unmapSharedMem((void *) n_ox, sizeof(*n_ox));
	unmapSharedMem((void *) barrier_cnt, sizeof(*barrier_cnt));
	sem_close(write_file);
	sem_close(h_q);
	sem_close(ox_q);
	sem_close(forming_h2o);
	sem_close(mutex);
	sem_close(reu_bar_turnstile1);
	sem_close(reu_bar_turnstile2);
	sem_unlink("/xjanot04.sem0");
	sem_unlink("/xjanot04.sem1");
	sem_unlink("/xjanot04.sem2");
	sem_unlink("/xjanot04.sem3");
	sem_unlink("/xjanot04.sem4");
	sem_unlink("/xjanot04.sem5");
	sem_unlink("/xjanot04.sem6");
}

int
init()
{
	/*
	 * initializing shared variables
	 */

	n_action = mapSharedMem(sizeof(*n_action));
	if (n_action == (void *) -1) {
		goto shm_err;
	}
	(*n_action) = 1;
	n_molecule = mapSharedMem(sizeof(*n_molecule));
	if (n_molecule == (void *) -1) {
		goto shm_err;
	}
	n_h = mapSharedMem(sizeof(*n_h));
	if (n_h == (void *) -1) {
		goto shm_err;
	}
	(*n_h) = 0;
	n_ox = mapSharedMem(sizeof(*n_ox));
	if (n_ox == (void *) -1) {
		goto shm_err;
	}
	(*n_ox) = 0;
	barrier_cnt = mapSharedMem(sizeof(*barrier_cnt));
	if (barrier_cnt == (void *) -1) {
		goto shm_err;
	}
	h_cnt = mapSharedMem(sizeof(*h_cnt));
	if (h_cnt == (void *) -1) {
		goto shm_err;
	}
	(*h_cnt) = num_h;
	ox_cnt = mapSharedMem(sizeof(*ox_cnt));
	if (ox_cnt == (void *) -1) {
		goto shm_err;
	}
	(*ox_cnt) = num_ox;

	/*
	 * initializing semaphores
	 */

	if ((write_file = sem_open("/xjanot04.sem0", O_CREAT | O_EXCL, 0666, 1)) == SEM_FAILED) {
		perror("sem_open");
		goto shm_err;
	}
	if ((h_q = sem_open("/xjanot04.sem1", O_CREAT | O_EXCL, 0666, 0)) == SEM_FAILED) {
		perror("sem_open");
		goto shm_err;
	}
	if ((ox_q = sem_open("/xjanot04.sem2", O_CREAT | O_EXCL, 0666, 0)) == SEM_FAILED) {
		perror("sem_open");
		goto shm_err;
	}
	if ((forming_h2o = sem_open("/xjanot04.sem3", O_CREAT | O_EXCL, 0666, 1)) == SEM_FAILED) {
		perror("sem_open");
		goto shm_err;
	}
	if ((mutex = sem_open("/xjanot04.sem4", O_CREAT | O_EXCL, 0666, 1)) == SEM_FAILED) {
		perror("sem_open");
		goto shm_err;
	}
	if ((reu_bar_turnstile1 = sem_open("/xjanot04.sem5", O_CREAT | O_EXCL, 0666, 0)) == SEM_FAILED) {
		perror("sem_open");
		goto shm_err;
	}
	if ((reu_bar_turnstile2 = sem_open("/xjanot04.sem6", O_CREAT | O_EXCL, 0666, 1)) == SEM_FAILED) {
		perror("sem_open");
		goto shm_err;
	}

	return 0;

	shm_err:
		cleanup();
		return -1;
}

int
main(int argc, char *argv[])
{
	srand(time(0));
	pid_t pid;

	if (parseArgs(argc, argv, &num_h, &num_ox, &max_queue_wait_time, &max_molecule_wait_time)) {
		return 1;
	}

	if (init()) {
		return 1;
	}

	/* created a file to write to if it does not already exist, if it does then truncate it */
	shared_file = fopen("proj2.out", "w+");
	if (!shared_file) {
		fprintf(stderr, "Couldn't open or create proj2.out\n");
		return 1;
	}

	/* loop over the amount of hydrogens and oxygens and fork them */
	for (int i = 1; i <= num_h; i++) {
		pid = fork();
		if (pid == 0) {
			hydrogen(i);
		}
		else if (pid == -1) {
			fprintf(stderr, "Error while forking hydrogen\n");
		}
	}
	for (int i = 1; i <= num_ox; i++) {
		pid = fork();
		if (pid == 0) {
			oxygen(i);
		}
		else if (pid == -1) {
			fprintf(stderr, "Error while forking oxygen\n");
		}
	}

	while (wait(NULL) > 0);
	cleanup();
	fclose(shared_file);
	return 0;
}