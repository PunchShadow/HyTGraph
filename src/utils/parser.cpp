/*
* Copyright 1997, Regents of the University of Minnesota
*
* Extracted from Metis io.c http://glaros.dtc.umn.edu/gkhome/metis/metis/download
*
* This file contains routines related to I/O
*
* Started 8/28/94
* George
*
* $Id: io.c 11932 2012-05-10 18:18:23Z dominique $
*
*/

#include <utils/parser.h>


//#include <crtdefs.h>
#include <iostream>
#include <fstream>
#include <sstream> 
#include <cstddef>
#include <cstdarg>
#include <cstring>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <set>
#include <vector>
#include <gflags/gflags.h>

static bool GetUndirectedFlag()
{
    std::string value;
    if (!gflags::GetCommandLineOption("undirected", &value))
        return false;

    if (value == "1" || value == "true" || value == "True" || value == "TRUE" ||
        value == "t" || value == "T" || value == "yes" || value == "Yes" || value == "YES")
        return true;

    return false;
}


struct Edge{
    idx_t src;
    idx_t dst;
};

struct EdgeWeighted{
    idx_t src;
    idx_t dst;
    idx_t weight;
};


/*************************************************************************/
/*! This function initializes a graph_t data structure */
/*************************************************************************/
void InitGraph(graph_t *graph)
{
    memset((void *) graph, 0, sizeof(graph_t));

    /* graph size constants */
    graph->nvtxs = -1;
    graph->nedges = -1;
    graph->ncon = -1;
    graph->mincut = -1;
    graph->minvol = -1;
    graph->nbnd = -1;

    /* memory for the graph structure */
    graph->xadj = NULL;
    graph->vwgt = NULL;
    graph->vsize = NULL;
    graph->adjncy = NULL;
    graph->adjwgt = NULL;
    graph->label = NULL;
    graph->cmap = NULL;
    graph->tvwgt = NULL;
    graph->invtvwgt = NULL;

    graph->readvw = false;
    graph->readew = false;

    /* by default these are set to true, but the can be explicitly changed afterwards */
    graph->free_xadj = 1;
    graph->free_vwgt = 1;
    graph->free_vsize = 1;
    graph->free_adjncy = 1;
    graph->free_adjwgt = 1;


    /* memory for the partition/refinement structure */
    graph->where = NULL;
    graph->pwgts = NULL;
    graph->id = NULL;
    graph->ed = NULL;
    graph->bndptr = NULL;
    graph->bndind = NULL;
    graph->nrinfo = NULL;
    graph->ckrinfo = NULL;
    graph->vkrinfo = NULL;

    /* linked-list structure */
    graph->coarser = NULL;
    graph->finer = NULL;
}

/*************************************************************************/
/*! This function creates and initializes a graph_t data structure */
/*************************************************************************/
graph_t *CreateGraph(void)
{
    graph_t *graph;

    graph = (graph_t *) malloc(sizeof(graph_t));

    InitGraph(graph);

    return graph;
}

/*************************************************************************/
/*! This function deallocates any memory stored in a graph */
/*************************************************************************/
void FreeGraph(graph_t **r_graph)
{
    graph_t *graph;

    graph = *r_graph;

    /* free graph structure */
    if (graph->free_xadj)
        free((void *) graph->xadj);
    if (graph->free_vwgt)
        free((void *) graph->vwgt);
    if (graph->free_vsize)
        free((void *) graph->vsize);
    if (graph->free_adjncy)
        free((void *) graph->adjncy);
    if (graph->free_adjwgt)
        free((void *) graph->adjwgt);

    /* free partition/refinement structure */
    //FreeRData(graph);

    free((void *) graph->tvwgt);
    free((void *) graph->invtvwgt);
    free((void *) graph->label);
    free((void *) graph->cmap);
    free((void *) graph);

    *r_graph = NULL;
}

//static int exit_on_error = 1;

/*************************************************************************/
/*! This function prints an error message and exits
*/
/*************************************************************************/
void errexit(const char *f_str, ...)
{
    va_list argp;

    va_start(argp, f_str);
    vfprintf(stderr, f_str, argp);
    va_end(argp);

    if (strlen(f_str) == 0 || f_str[strlen(f_str) - 1] != '\n')
        fprintf(stderr, "\n");
    fflush(stderr);

    if (/*exit_on_error*/ 1)
        exit(-2);

    /* abort(); */
}

/*************************************************************************
* This function opens a file
**************************************************************************/
FILE *gk_fopen(const char *fname, const char *mode, const char *msg)
{
    FILE *fp;
    char errmsg[8192];

    fp = fopen(fname, mode);
    if (fp != NULL)
        return fp;

    sprintf(errmsg, "file: %s, mode: %s, [%s]", fname, mode, msg);
    perror(errmsg);
    errexit("Failed on gk_fopen()\n");

    return NULL;
}


/*************************************************************************
* This function closes a file
**************************************************************************/
void gk_fclose(FILE *fp)
{
    fclose(fp);
}


/*************************************************************************/
/*! This function is the GKlib implementation of glibc's getline()
function.
\returns -1 if the EOF has been reached, otherwise it returns the
number of bytes read.
*/
/*************************************************************************/
ptrdiff_t gk_getline(char **lineptr, size_t *n, FILE *stream)
{
#ifdef HAVE_GETLINE
    return getline(lineptr, n, stream);
#else
    size_t i;
    int ch;

    if (feof(stream))
        return -1;

    /* Initial memory allocation if *lineptr is NULL */
    if (*lineptr == NULL || *n == 0)
    {
        *n = 1024;
        *lineptr = (char *) malloc((*n) * sizeof(char));
    }

    /* get into the main loop */
    i = 0;
    while ((ch = getc(stream)) != EOF)
    {
        (*lineptr)[i++] = (char) ch;

        /* reallocate memory if reached at the end of the buffer. The +1 is for '\0' */
        if (i + 1 == *n)
        {
            *n = 2 * (*n);
            *lineptr = (char *) realloc(*lineptr, (*n) * sizeof(char));
        }

        if (ch == '\n')
            break;
    }
    (*lineptr)[i] = '\0';

    return (i == 0 ? -1 : i);
#endif
}

/*************************************************************************/
/*! This function reads in a sparse graph */
/*************************************************************************/
graph_t *ReadGraph(char *filename)
{
    idx_t i, j, k, l, fmt, ncon, nfields, readew, readvw, readvs, edge, ewgt;
    idx_t *xadj, *adjncy, *vwgt, *adjwgt, *vsize;
    char *line = NULL, fmtstr[256], *curstr, *newstr;
    size_t lnlen = 0;
    FILE *fpin;
    graph_t *graph;

    graph = CreateGraph();


    return graph;
}

#ifdef WIN32
// Windows "host" byte order is little endian
static inline uint64_t le64toh(uint64_t x) {
    return x;
}

#endif

/*************************************************************************/
/*! This function reads in a sparse graph */
/*************************************************************************/
graph_t *ReadGraphGR(char *filename)
{
    idx_t *xadj, *adjncy, *vwgt, *adjwgt, *vsize;
    FILE *fpin;
    graph_t *graph;

    graph = CreateGraph();
    return graph;
}

graph_t *ReadGraphMarket(char *filename)
{
    idx_t *xadj, *adjncy, *vwgt, *adjwgt, *vsize;
    FILE *fpin;
    char *line = NULL;
    size_t lnlen = 0;
    graph_t *graph;

    graph = CreateGraph();


    return graph;
}


graph_t *ReadGraphMarket_bigdata(char *filename,idx_t weight_num)
{
    uint64_t *xadj;
    idx_t *adjncy, *vwgt, *adjwgt, *vsize;

    std::ifstream infile;
    infile.open(filename);
    std::stringstream ss;
    std::string line;

    graph_t *graph;

    graph = CreateGraph();

    int weighted = 1; // 0: no weight; 1: int weight; 2: generated weight

    if(weight_num == 1){
        weighted = 2;
    }

    graph->nedges = 0;
    graph->nvtxs = 0;

    std::vector<uint64_t> xadj_pri;
    
    xadj_pri.resize(1);

    uint32_t src = 0, dst = 0;
    const bool undirected_flag = GetUndirectedFlag();
    uint64_t max_node = 0;
    uint64_t header_nodes = 0;
    uint64_t header_edges = 0;
    bool header_found = false;
    bool header_undirected = false;

    auto parse_header_counts = [&](const std::string &l) {
        unsigned long long n = 0, e = 0;
        if (std::sscanf(l.c_str(), "# Nodes: %llu Edges: %llu", &n, &e) == 2)
        {
            header_nodes = n;
            header_edges = e;
            header_found = true;
        }
        if (l.find("Undirected") != std::string::npos)
        {
            header_undirected = true;
        }
    };

    while(getline( infile, line )){
        const size_t pos = line.find_first_not_of(" \t\r\n");
        if (pos == std::string::npos)
            continue;
        if (line[pos] == '#')
        {
            parse_header_counts(line);
            continue;
        }

        ss.str("");
        ss.clear();
        ss << line;

        if (!(ss >> src >> dst))
            continue;

        const uint64_t local_max = (src > dst) ? src : dst;
        if (max_node < local_max)
            max_node = local_max;

        graph->nedges++;

        if (xadj_pri.size() <= src)
        {
            size_t new_size = xadj_pri.size() > 0 ? xadj_pri.size() : 1;
            while (new_size <= src)
                new_size *= 2;
            xadj_pri.resize(new_size, 0);
        }
        xadj_pri[src]++;

        if (undirected_flag && src != dst)
        {
            if (xadj_pri.size() <= dst)
            {
                size_t new_size = xadj_pri.size() > 0 ? xadj_pri.size() : 1;
                while (new_size <= dst)
                    new_size *= 2;
                xadj_pri.resize(new_size, 0);
            }
            xadj_pri[dst]++;
            graph->nedges++;
        }
    }
    infile.close();

    uint64_t nvtxs = max_node + 1;
    if (header_found && header_nodes > nvtxs)
        nvtxs = header_nodes;
    graph->nvtxs = static_cast<idx_t>(nvtxs);

    if (xadj_pri.size() < graph->nvtxs)
    {
        xadj_pri.resize(graph->nvtxs, 0);
    }

    if (header_undirected && !undirected_flag)
    {
        fprintf(stderr, "Note: SNAP header indicates undirected graph. "
                        "Use --undirected to symmetrize edges.\n");
    }

    //vwgt = graph->vwgt = (idx_t *) calloc((0 * graph->nvtxs), sizeof(idx_t));  // file doesn't store node weights though.
    graph->readvw = false;

    xadj = graph->xadj = (uint64_t *) calloc((graph->nvtxs + 1), sizeof(uint64_t));
    if(weighted == 0){
        adjncy = graph->adjncy = (idx_t *) calloc((graph->nedges), sizeof(uint32_t));
    }
    else{
        graph->readew = true;
        adjncy = graph->adjncy = (idx_t *) calloc((graph->nedges), sizeof(uint32_t));
        adjwgt = graph->adjwgt = (idx_t *) calloc((graph->nedges), sizeof(uint32_t));
    }

    uint64_t count = 0;
    for (idx_t src = 0; src < graph->nvtxs; src++)
    {
        xadj[src] = count;
        count += xadj_pri[src];
    }
    xadj[graph->nvtxs] = graph->nedges;
    

    infile.open(filename);

    idx_t *outDegreeCounter  = (idx_t *) calloc((graph->nvtxs + 1), sizeof(idx_t));
    for(idx_t i=0; i<graph->nvtxs; i++)
        outDegreeCounter[i] = 0;
    uint32_t weight;
    bool warned_missing_weight = false;
    while(getline( infile, line )){
        
        const size_t pos = line.find_first_not_of(" \t\r\n");
        if (pos == std::string::npos)
            continue;
        if (line[pos] == '#')
            continue;

        ss.str("");
        ss.clear();
        ss << line;

        if (!(ss >> src >> dst))
            continue;

 
        uint64_t location = xadj[src] + outDegreeCounter[src];                
        adjncy[location] = dst;
        outDegreeCounter[src]++; 

        uint32_t edge_weight = 1;

        if(weighted == 2){
            adjwgt[location] = src % 64;
        }
        else if(weighted == 1){
            if (ss >> weight)
            {
                edge_weight = weight;
                adjwgt[location] = edge_weight;
            }
            else
            {
                if (!warned_missing_weight)
                {
                    fprintf(stderr, "Warning: missing edge weights; defaulting to 1\n");
                    warned_missing_weight = true;
                }
                adjwgt[location] = edge_weight;
            }
        }

        if (undirected_flag && src != dst)
        {
            uint64_t rev_location = xadj[dst] + outDegreeCounter[dst];
            adjncy[rev_location] = src;
            outDegreeCounter[dst]++;

            if (weighted == 2)
            {
                adjwgt[rev_location] = dst % 64;
            }
            else if (weighted == 1)
            {
                adjwgt[rev_location] = edge_weight;
            }
        }
    }
    infile.close();
    free(outDegreeCounter);

    return graph;
}

graph_t *ReadGraphBCSR(char *filename, bool weighted)
{
    std::ifstream infile(filename, std::ios::in | std::ios::binary);
    if (!infile.is_open()) {
        errexit("Failed to open binary CSR file: %s\n", filename);
    }

    uint32_t num_nodes = 0;
    uint32_t num_edges = 0;
    infile.read(reinterpret_cast<char*>(&num_nodes), sizeof(uint32_t));
    infile.read(reinterpret_cast<char*>(&num_edges), sizeof(uint32_t));
    if (!infile) {
        errexit("Failed to read binary CSR header: %s\n", filename);
    }

    graph_t *graph = CreateGraph();
    graph->nvtxs = static_cast<idx_t>(num_nodes);
    graph->nedges = num_edges;
    graph->readvw = false;
    graph->readew = weighted;

    graph->xadj = (uint64_t *) calloc((graph->nvtxs + 1), sizeof(uint64_t));
    graph->adjncy = (idx_t *) calloc((graph->nedges), sizeof(uint32_t));
    if (weighted) {
        graph->adjwgt = (idx_t *) calloc((graph->nedges), sizeof(uint32_t));
    }

    std::vector<uint32_t> row_offsets(num_nodes);
    if (num_nodes > 0) {
        infile.read(reinterpret_cast<char*>(row_offsets.data()), sizeof(uint32_t) * num_nodes);
        if (!infile) {
            errexit("Failed to read binary CSR row offsets: %s\n", filename);
        }
    }
    for (uint32_t i = 0; i < num_nodes; i++) {
        graph->xadj[i] = row_offsets[i];
    }
    graph->xadj[num_nodes] = num_edges;

    if (weighted) {
        struct BcsrEdgeWeighted {
            uint32_t end;
            uint32_t w8;
        };
        std::vector<BcsrEdgeWeighted> edges(num_edges);
        if (num_edges > 0) {
            infile.read(reinterpret_cast<char*>(edges.data()), sizeof(BcsrEdgeWeighted) * num_edges);
            if (!infile) {
                errexit("Failed to read binary CSR weighted edges: %s\n", filename);
            }
        }
        for (uint32_t i = 0; i < num_edges; i++) {
            graph->adjncy[i] = static_cast<idx_t>(edges[i].end);
            graph->adjwgt[i] = static_cast<idx_t>(edges[i].w8);
        }
    } else {
        if (num_edges > 0) {
            infile.read(reinterpret_cast<char*>(graph->adjncy), sizeof(uint32_t) * num_edges);
            if (!infile) {
                errexit("Failed to read binary CSR edges: %s\n", filename);
            }
        }
    }

    infile.close();
    return graph;
}
