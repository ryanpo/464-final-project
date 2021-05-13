using System.Collections;
using System.Collections.Generic;
using UnityEngine;


public class FEM_cube : MonoBehaviour
{
    public class Node
    {
        public Vector3 position;
        public Vector3 velocity;
        public bool is_stationary;
        
        public Node(Vector3 pos, Vector3 vel, bool stationary = false) {
            position = pos;
            is_stationary = false;
            velocity = vel;
        }

    }

    public class Spring
    {
        public int node1;
        public int node2;
        public float k;
        public float x;
        public bool active;
        public Vector3 force;

        public Spring(int n1, int n2, float sc, float rest_len, bool act = true) {
            node1 = n1;
            node2 = n2;
            k = sc;
            x = rest_len;
            active = act;
        }
    }


    public float lamb = 11.1f;
    public float mu = 11.1f;
    public float phi = 0.00f;
    public float psi = 0.00f;
    public float dt = 0.1f;
    public Transform sphere_trans;

    [HideInInspector]
    public Vector3[] vertices_m;
    [HideInInspector]
    public Vector3[] vertices_p;
    [HideInInspector]
    public int[] tris;
    [HideInInspector]
    public Vector3[] normals;
    [HideInInspector]
    public Vector2[] uv;
    [HideInInspector]
    public Node[] nodes;
    [HideInInspector]
    public Spring[] springs;
    [HideInInspector]
    public MeshRenderer meshRenderer;
    [HideInInspector]
    public MeshFilter meshFilter;
    [HideInInspector]
    public Mesh mesh;
    [HideInInspector]
    public Vector3[] forces;
    [HideInInspector]
    public int net_count = 0;
	[HideInInspector]
	public Color[] colors;
	[HideInInspector]
	public Material material;
	[HideInInspector]
	public Vector3[] velocities;
	[HideInInspector]
	public int[][] tets;

    public void Init() {

    	material = GetComponent<Renderer> ().material;
    	material.SetColor("_Color", Color.blue);
    	material.SetFloat("_Metallic", 0.3f);
 		material.SetFloat("_GlossMapScale", 0.43f);
    	dt = 0.1f;
    	lamb = 5.0f;
    	mu = 5.0f;
    	phi = 3.0f;
    	psi = 3.0f;
    	vertices_m = new Vector3[8];
    	vertices_m[0] = new Vector3(1.0f, 0.0f, 1.0f);
    	vertices_m[1] = new Vector3(-1.0f, 0.0f, 1.0f);
    	vertices_m[2] = new Vector3(-1.0f, 0.0f, -1.0f);
    	vertices_m[3] = new Vector3(1.0f, 0.0f, -1.0f);	
    	vertices_m[4] = new Vector3(1.0f, 2.0f, 1.0f);
    	vertices_m[5] = new Vector3(-1.0f, 2.0f, 1.0f);
    	vertices_m[6] = new Vector3(-1.0f, 2.0f, -1.0f);
    	vertices_m[7] = new Vector3(1.0f, 2.0f, -1.0f);	

    	vertices_p = new Vector3[8];
    	vertices_p[0] = new Vector3(1.0f, 3.0f, 1.0f);
    	vertices_p[1] = new Vector3(-1.0f, 3.0f, 1.0f);
    	vertices_p[2] = new Vector3(-1.0f, 3.0f, -1.0f);
    	vertices_p[3] = new Vector3(1.0f, 3.0f, -1.0f);	
    	vertices_p[4] = new Vector3(1.0f, 5.0f, 1.0f);
    	vertices_p[5] = new Vector3(-1.0f, 5.0f, 1.0f);
    	vertices_p[6] = new Vector3(-1.0f, 5.0f, -1.0f);
    	vertices_p[7] = new Vector3(1.0f, 5.0f, -1.0f);	

    	for (int i = 0; i < 8; i++) {
    		vertices_p[i] += new Vector3(0.0f, 5.0f, 0.0f);
    	}

    	forces = new Vector3[8];
    	velocities = new Vector3[8];
    	velocities[0] = new Vector3(0.0f, 0.0f, 0.0f);
    	velocities[1] = new Vector3(0.0f, 0.0f, 0.0f);
    	velocities[2] = new Vector3(0.0f, 0.0f, 0.0f);
    	velocities[3] = new Vector3(0.0f, 0.0f, 0.0f);
    	velocities[4] = new Vector3(0.0f, 0.0f, 0.0f);
    	velocities[5] = new Vector3(0.0f, 0.0f, 0.0f);
    	velocities[6] = new Vector3(0.0f, 0.0f, 0.0f);
    	velocities[7] = new Vector3(0.0f, 0.0f, 0.0f);

    	tets = new int[5][];
    	tets[0] = new int[4];
    	tets[0][0] = 0;
    	tets[0][1] = 2;
    	tets[0][2] = 3;
    	tets[0][3] = 7;
    	
    	tets[1] = new int[4];
    	tets[1][0] = 0;
    	tets[1][1] = 2;
    	tets[1][2] = 7;
    	tets[1][3] = 5;
    	
    	tets[2] = new int[4];
    	tets[2][0] = 0;
    	tets[2][1] = 5;
    	tets[2][2] = 7;
    	tets[2][3] = 4;
    	
    	tets[3] = new int[4];
    	tets[3][0] = 2;
    	tets[3][1] = 7;
    	tets[3][2] = 5;
    	tets[3][3] = 6;

    	tets[4] = new int[4];
    	tets[4][0] = 0;
    	tets[4][1] = 1;
    	tets[4][2] = 2;
    	tets[4][3] = 5;

    	tris = new int[5*4*3];

    	for (int i = 0; i < 5; i++) {
    		tris[12*i    ] = tets[i][0];
    		tris[12*i + 1] = tets[i][1];
    		tris[12*i + 2] = tets[i][2];

    		tris[12*i + 3] = tets[i][3];
    		tris[12*i + 4] = tets[i][2];
    		tris[12*i + 5] = tets[i][1];

    		tris[12*i + 6] = tets[i][0];
    		tris[12*i + 7] = tets[i][3];
    		tris[12*i + 8] = tets[i][1];

    		tris[12*i + 9] = tets[i][0];
    		tris[12*i + 10] =tets[i][2];
    		tris[12*i + 11] =tets[i][3];
    	}


    	normals = new Vector3[8];
    	for (int i = 0; i < 8; i++) {
    		normals[i] = -Vector3.forward;
    	}
    	uv = new Vector2[8];
    	material = GetComponent<Renderer> ().material;
    	material.SetColor("_Color", Color.blue);
    	material.SetFloat("_Metallic", 0.6f);
 		material.SetFloat("_GlossMapScale", 0.43f);
    	mesh.vertices = vertices_p;
        mesh.triangles = tris;
        mesh.normals = normals;
        mesh.uv = uv;
        meshFilter.mesh = mesh;

    }

    private void DrawLines(int[] tetv) {
    	for (int i = 0; i < 4; i++) {
    		for (int j = 0; j < 4; j++) {
    			if (j < i) {
        			Debug.DrawLine(vertices_p[tetv[i]], vertices_p[tetv[j]], Color.blue, Time.deltaTime);
    			}
    		}
    	}
    }


    private void FixedUpdate() {
    	ClearForces();
    	for (int i = 0; i < 5; i++) {
    		GetForces(tets[i]);
    	}
    	GetPositions();
    	mesh.vertices = vertices_p;
    	for (int i = 0; i < 5; i++) {
    		DrawLines(tets[i]);
    	}
    	
        return;
    }

    private void ClearForces() {
    	for (int i = 0; i < 8; i++) {
    		forces[i] = new Vector3(0.0f, 0.0f, 0.0f);
    	}
    }

    private void GetForces(int[] tetvs) {
    	Matrix4x4 beta_inv = new Matrix4x4();
    	for (int i= 0; i < 4; i++) {
    		beta_inv.SetColumn(i, new Vector4(vertices_m[tetvs[i]].x, vertices_m[tetvs[i]].y, vertices_m[tetvs[i]].z, 1.0f));
    	}

    	Matrix4x4 beta = beta_inv.inverse;


    	Vector3[] dx_du = new Vector3[3];
    	for (int i = 0; i < 3; i++) {
    		dx_du[i] = dxdu_helper(beta, i, tetvs);
    	}

    	Vector3[] dxp_du = new Vector3[3];
    	for (int i = 0; i < 3; i++) {
    		dxp_du[i] = dxpdu_helper(beta, i, tetvs);
    	}

    	Matrix4x4 strain = new Matrix4x4();
    	for (int i = 0; i < 3; i++) {
    		for (int j = 0; j < 3; j++) {
    			if (i == j) {
    				strain[i, j] = Vector3.Dot(dx_du[i], dx_du[j]) - 1.0f;
    			} else {
    				strain[i, j] = Vector3.Dot(dx_du[i], dx_du[j]);
    			}
    		}
    	}

    	Matrix4x4 strain_rate = new Matrix4x4();
    	for (int i = 0; i < 3; i++) {
    		for (int j = 0; j < 3; j++) {
    			strain_rate[i, j] = Vector3.Dot(dx_du[i], dxp_du[j]) + Vector3.Dot(dx_du[j], dxp_du[i]);
    		}
    	}


    	Matrix4x4 sig_e = new Matrix4x4();
    	for (int i = 0; i < 3; i++) {
    		for (int j = 0; j < 3; j++) {
    			sig_e[i, j] = 0;
    			for (int k = 0; k < 3; k++) {
    				if (i == j) {
    					sig_e[i, j] += 2*mu*strain[i, j];
    				} else {
    					sig_e[i, j] += lamb*strain[k, k] + 2*mu*strain[i, j];
    				}
    			}
    		}
    	}

    	

    	for (int i = 0; i < 3; i++) {
    		for (int j = 0; j < 3; j++) {
    			for (int k = 0; k < 3; k++) {
    				if (i == j) {
    					sig_e[i, j] += 2*psi*strain_rate[i, j];
    				} else {
    					sig_e[i, j] += phi*strain_rate[k, k] + 2*psi*strain_rate[i, j];
    				}
    			}
    		}
    	}



    	float vol = Vector3.Dot(Vector3.Cross(vertices_m[tetvs[1]] - vertices_m[tetvs[0]], 
    										  vertices_m[tetvs[2]] - vertices_m[tetvs[0]]),
    										  vertices_m[tetvs[3]] - vertices_m[tetvs[0]]);
    	vol = vol/6.0f;

    	for (int i = 0; i < 4; i++) {
    		for (int j = 0; j < 4; j++) {
    			float tmp = 0.0f;
    			for (int k = 0; k < 3; k++) {
    				for (int l = 0; l < 3; l++) {
    					tmp += beta[j, l] * beta[i, k] * sig_e[k, l];
    				}
    			}
    			forces[tetvs[i]] += tmp*vertices_p[tetvs[j]];
    		}
    		forces[tetvs[i]] = forces[tetvs[i]] * (-vol/2.0f);
    	}
    }

    private Vector3 dxdu_helper(Matrix4x4 beta, int i, int[] tetvs) {
    	Vector3 ret = new Vector3();
    	for (int j = 0; j < 3; j++) {
    		ret[j] = Vector4.Dot(beta.GetColumn(i), new Vector4(vertices_p[tetvs[0]][j], vertices_p[tetvs[1]][j], 
    														vertices_p[tetvs[2]][j], vertices_p[tetvs[3]][j]));
    	}
    	return ret;
    }

    private Vector3 dxpdu_helper(Matrix4x4 beta, int i, int[] tetvs) {
    	Vector3 ret = new Vector3();
    	for (int j = 0; j < 3; j++) {
    		ret[j] = Vector4.Dot(beta.GetColumn(i), new Vector4(velocities[tetvs[0]][j], velocities[tetvs[1]][j], 
    														velocities[tetvs[2]][j], velocities[tetvs[3]][j]));
    	}
    	return ret;
    }

    public void Start() {
    	meshRenderer = gameObject.AddComponent<MeshRenderer>();
        meshRenderer.sharedMaterial = new Material(Shader.Find("Standard"));
        meshFilter = gameObject.AddComponent<MeshFilter>();
        mesh = new Mesh();

    	Init();
    }

    private void GetPositions() {

    	for (int i = 0; i < 8; i++) {
    		forces[i] += new Vector3(0.0f, 0.2f, 0.0f);
    	}

    	for (int i = 0; i < 8; i++) {
    		float sphere_r = sphere_trans.GetComponent<SphereCollider>().radius * sphere_trans.transform.localScale.x;
	        Vector3 sphere_pos = sphere_trans.position;
	        float cc = 2.4f;
	        float cr = sphere_r;

	        if ((vertices_p[i] - sphere_pos).magnitude < cr) {
                    
                Vector3 diff = vertices_m[i] - sphere_pos;
                forces[i] -= (cr-diff.magnitude) * (diff/diff.magnitude) * cc;
            }

    		velocities[i] -= forces[i]/30;
    		vertices_p[i] += velocities[i]*dt;
    	}



    	//Debug.Log(velocities[0].ToString("F10"));
    	for (int i = 0; i < 8; i++) {
    		if (vertices_p[i].y < 0) {
    			vertices_p[i].y = 0.0f;
    			velocities[i] = new Vector3(0.0f,0.0f,0.0f);
    		}
    	}
    	return;
    }
    
}


