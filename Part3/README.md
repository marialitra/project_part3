# project_part3

1.)

ALSO CHANGE THE FOLLOWING IN PART 1 POF THE PROJECT WHICH WE WILL UPLOAD IN GITHUB!!!!!
IMPORTANT CHANGE!!!!!!!!!!!!!!!

changes in ivfpq c source code, in ivfpq_init about 
    // index->k = k_clusters; line 335

    and change it into:
    
    if(ivf_temp->k <= 0)
    {
        fprintf(stderr, "Error in kclusters\n");
        exit(1);
    }
    
    index->k = ivf_temp->k;
    k_clusters = index->k;

    lines 351 - 358


2.)

not called ids.txt but verctors_ids.txt, problem?

3.) in query.c of the first exercise we re entered the lines about the count of the QPS
    and we write it at the end of the file so that python part3 handles it as it wants to

4.) in nlsh_core.py from part2 we inserted the things about computing the QPS (there were there before)

5.) problem with M in pq and hypercube share the same in defaults so we shall see

6.) BLAST takes the same N as search will take