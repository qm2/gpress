#ifndef HASH_H
#define HASH_H



//Hashtable element structure
typedef struct hash_elem_t {
	struct hash_elem_t* next; // Next element in case of a collision
	void* data;	// Pointer to the stored element
	char key[]; 	// Key of the stored element
} hash_elem_t;

//Hashtabe structure
typedef struct {
	unsigned int capacity;	// Hashtable capacity (in terms of hashed keys)
	unsigned int e_num;	// Number of element currently stored in the hashtable
	hash_elem_t** table;	// The table containaing elements
} hashtable_t;

//Structure used for iterations
typedef struct {
	hashtable_t* ht; 	// The hashtable on which we iterate
	unsigned int index;	// Current index in the table
	hash_elem_t* elem; 	// Curent element in the list
} hash_elem_it;


static unsigned int hash_table_calc_hash(char* str);

hashtable_t* ht_create(unsigned int capacity);

void* ht_put(hashtable_t* hasht, char* key, void* data);

void* ht_get(hashtable_t* hasht, char* key);

void* ht_remove(hashtable_t* hasht, char* key);

void ht_list_keys(hashtable_t* hasht, char** k, size_t len);

void ht_list_values(hashtable_t* hasht, void** v, size_t len);

hash_elem_t* ht_iterate(hash_elem_it* iterator);

char* ht_iterate_keys(hash_elem_it* iterator);

void* ht_iterate_values(hash_elem_it* iterator);

void ht_clear(hashtable_t* hasht, int free_data);

void ht_destroy(hashtable_t* hasht);

#endif
