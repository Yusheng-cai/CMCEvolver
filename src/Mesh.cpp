#include "Mesh.h"

Mesh::Mesh(const std::vector<Real3>& v, const std::vector<INT3>& f){
    vertices_.clear();
    triangles_.clear();

    for (int i=0;i<v.size();i++){
        vertex vv;
        vv.position_ = v[i];
        vertices_.push_back(vv);
    }

    triangles_.clear();
    for (int i=0;i<f.size();i++){
        triangle t;
        t.triangleindices_ = f[i];
        triangles_.push_back(t);
    }

    CalcVertexNormals();
}

void Mesh::MoveVertexIntoBox(const Real3& OldVerPos, Real3& NewVerPos)
{
    if (isPeriodic()){
        Real3 boxCenter;
        boxCenter = boxLength_ * 0.5;

        Real3 diff;
        for (int i=0;i<3;i++){
            diff[i] = OldVerPos[i] - boxCenter[i];

            if (diff[i] >= boxLength_[i] * 0.5) diff[i] = diff[i] - boxLength_[i];
            if (diff[i] <= -boxLength_[i] * 0.5) diff[i] = diff[i] + boxLength_[i]; 

            NewVerPos[i] = boxCenter[i] + diff[i];
        }
    }
}

void Mesh::SetVerticesAndTriangles(const std::vector<vertex>& v, const std::vector<triangle>& t){
    vertices_.clear();
    vertices_.insert(vertices_.end(),v.begin(), v.end());

    triangles_.clear();
    triangles_.insert(triangles_.end(), t.begin(), t.end());
}

Mesh::Real3 Mesh::CalculateCOM(){
    Real3 COM= {{0,0,0}};
    for (auto v : vertices_){
        COM = COM + v.position_;
    }

    COM = COM / vertices_.size();

    return COM;
}

void Mesh::ShiftCOMWithRespectTo(Real3& COM)
{
    Real3 COM_current = CalculateCOM();
    for (auto& v : vertices_){
        v.position_ = v.position_ - COM_current + COM;
    }
}

void Mesh::CalcPerVertexDir()
{
    PerVertexdir1_.resize(vertices_.size());
    PerVertexdir2_.resize(vertices_.size());

    for (int i=0;i<triangles_.size();i++){
        auto& t = triangles_[i];
 
        int index0 = t[0];
        int index1 = t[1];
        int index2 = t[2];

        Real3 diff0, diff1, diff2;
        Real diffsq0, diffsq1, diffsq2;

        getVertexDistance(vertices_[index1], vertices_[index0], diff0, diffsq0);
        getVertexDistance(vertices_[index2], vertices_[index1], diff1, diffsq1);
        getVertexDistance(vertices_[index0], vertices_[index2], diff2, diffsq2);

        PerVertexdir1_[index0] = diff0;
        PerVertexdir1_[index1] = diff1;
        PerVertexdir1_[index2] = diff2;
    }

    for (int i=0;i<PerVertexdir1_.size();i++)
    {
        PerVertexdir1_[i] = LinAlg3x3::CrossProduct(PerVertexdir1_[i], vertices_[i].normals_);
        LinAlg3x3::normalize(PerVertexdir1_[i]);

        Real3 B = LinAlg3x3::CrossProduct(vertices_[i].normals_, PerVertexdir1_[i]) ;
        LinAlg3x3::normalize(B);

        PerVertexdir2_[i] = B;
    } 
}

Mesh::Real Mesh::calculateVolume()
{
    Real volume_ =0.0;
    for (int i=0;i<triangles_.size();i++)
    {
        auto& t = triangles_[i];

        Real3 g;
        g.fill(0);
        for (int j=0;j<3;j++)
        {
            int index= t[j];
            for (int k=0;k<3;k++){
                g[k] += vertices_[index].position_[k] /3.0;
            }
        }

        Real3 vec1, vec2;
        int index1 = t[0];
        int index2 = t[1];
        int index3 = t[2];

        for (int j=0;j<3;j++){
            vec1[j] = vertices_[index2].position_[j] - vertices_[index1].position_[j];
            vec2[j] = vertices_[index3].position_[j] - vertices_[index1].position_[j];
        }

        Real3 N = LinAlg3x3::CrossProduct(vec1, vec2);

        volume_ += 1.0/6.0 * LinAlg3x3::DotProduct(g,N);
    }

    return volume_;
}

void Mesh::scaleVertices(Real num)
{
    #pragma omp parallel for
    for (int i=0;i<vertices_.size();i++){
        vertices_[i].position_ = num * vertices_[i].position_;
    }

    update();
}

MeshTools::Real MeshTools::CalculateVolumeUnderneath(Mesh& m, int projected_plane){
    // obtain the vertices and triangles 
    const auto& vertices = m.getvertices();
    const auto& faces    = m.gettriangles();

    std::vector<Real> vecAreas;
    std::vector<Real3> Normals;
    MeshTools::CalculateTriangleAreasAndFaceNormals(m, vecAreas, Normals);

    Real volume = 0.0;

    #pragma omp parallel
    {
        Real v_local = 0.0;

        #pragma omp for
        for (int i=0;i<faces.size();i++){
            // get the vertices
            auto v1 = vertices[faces[i].triangleindices_[0]];
            auto v2 = vertices[faces[i].triangleindices_[1]];
            auto v3 = vertices[faces[i].triangleindices_[2]];

            // average height
            Real avg_height = 1.0 / 3.0 * (v1.position_[projected_plane] + v2.position_[projected_plane] + v3.position_[projected_plane] 
                                            );
            v_local += vecAreas[i] * Normals[i][projected_plane] * avg_height;
        }

        #pragma omp critical
        {
            volume += v_local;
        }
    }

    return volume;
}


void Mesh::getVertexDistance(const Real3& v1, const Real3& v2, Real3& distVec, Real& dist) const 
{
    distVec.fill(0);
    distVec = v1 - v2;

    if (isPeriodic()){
        for (int i=0;i<3;i++){
            if (distVec[i] > boxLength_[i] * 0.5) distVec[i] -= boxLength_[i];
            if (distVec[i] < -boxLength_[i] * 0.5) distVec[i] += boxLength_[i];
        }
    }
    dist = LinAlg3x3::norm(distVec);
}

void Mesh::getVertexDistance(const vertex& v1, const vertex& v2, Real3& distVec, Real& dist) const
{
    getVertexDistance(v1.position_, v2.position_, distVec, dist);
}

void Mesh::CalculateShift(const Real3& v1, const Real3& v2, Real3& shiftVec)
{
    for (int i=0;i<3;i++){
        Real d = v1[i] - v2[i];

        if (d > (0.5 * boxLength_[i])) shiftVec[i] = -boxLength_[i];
        else if (d < (- 0.5 * boxLength_[i])) shiftVec[i] = boxLength_[i];
        else shiftVec[i] = 0.0;
    }
}

Mesh::Real3 Mesh::getShiftedVertexPosition(const vertex& v1, const vertex& v2) const
{
    Real3 dist;
    Real distsq;

    // v1 - v2
    getVertexDistance(v1, v2, dist, distsq);
    Real3 ret = v2.position_ + dist;

    return ret;
}

Mesh::Real3 Mesh::getShiftedVertexPosition(const Real3& v1, const Real3& v2) const
{
    Real3 dist;
    Real distsq;

    // v1 - v2
    getVertexDistance(v1, v2, dist, distsq);
    Real3 ret = v2 + dist;

    return ret;
}

Mesh::Real3 Mesh::getShiftIntoBox(const Real3& v1)
{
    Real3 shift;
    Real3 center;
    center = boxLength_ * 0.5;

    for (int i=0;i<3;i++){
        Real diff = v1[i] - center[i];

        if (diff >= boxLength_[i] * 0.5) shift[i] = -boxLength_[i];
        else if (diff <= - boxLength_[i] * 0.5) shift[i] = boxLength_[i];
        else shift[i] = 0.0;
    }

    return shift;
}

// Compute per-vertex point areas
void Mesh::CalculateCornerArea()
{
    triangleArea_.clear();
    triangleArea_.resize(vertices_.size(),0.0);

    int nf = triangles_.size();

    cornerArea_.clear();
    cornerArea_.resize(nf);

	for (int i = 0; i < nf; i++) {
		// Edges
        auto& t = triangles_[i];
        int index1 = t[0];
        int index2 = t[1];
        int index3 = t[2];

        Real3 edge1, edge2, edge3;
        Real edge1sq, edge2sq, edge3sq;

        getVertexDistance(vertices_[index3], vertices_[index2], edge1, edge1sq);
        getVertexDistance(vertices_[index1], vertices_[index3], edge2, edge2sq);
        getVertexDistance(vertices_[index2], vertices_[index1], edge3, edge3sq);

		// Compute corner weights
        Real3 cross = LinAlg3x3::CrossProduct(edge1, edge2);
        Real  area  = 0.5 * LinAlg3x3::norm(cross);

        Real e1 = LinAlg3x3::norm(edge1);
        Real e2 = LinAlg3x3::norm(edge2);
        Real e3 = LinAlg3x3::norm(edge3);

		Real3 l2 = { e1*e1, e2*e2, e3*e3};

		// Barycentric weights of circumcenter
		Real3 bcw = { l2[0] * (l2[1] + l2[2] - l2[0]),
		                 l2[1] * (l2[2] + l2[0] - l2[1]),
		                 l2[2] * (l2[0] + l2[1] - l2[2]) };

		if (bcw[0] <= 0.0f) {
			cornerArea_[i][1] = -0.25f * l2[2] * area/LinAlg3x3::DotProduct(edge1, edge3);
			cornerArea_[i][2] = -0.25f * l2[1] * area/LinAlg3x3::DotProduct(edge1, edge2);
			cornerArea_[i][0] = area - cornerArea_[i][1] - cornerArea_[i][2];
		} else if (bcw[1] <= 0.0f) {
			cornerArea_[i][2] = -0.25f * l2[0] * area/LinAlg3x3::DotProduct(edge2, edge1);
			cornerArea_[i][0] = -0.25f * l2[2] * area/LinAlg3x3::DotProduct(edge2, edge3);
			cornerArea_[i][1] = area - cornerArea_[i][2] - cornerArea_[i][0];
		} else if (bcw[2] <= 0.0f) {
			cornerArea_[i][0] = -0.25f * l2[1] * area/LinAlg3x3::DotProduct(edge3, edge2);
			cornerArea_[i][1] = -0.25f * l2[0] * area/LinAlg3x3::DotProduct(edge3, edge1);
			cornerArea_[i][2] = area - cornerArea_[i][0] - cornerArea_[i][1];
		} else {
			float scale = 0.5f * area / (bcw[0] + bcw[1] + bcw[2]);
			for (int j = 0; j < 3; j++)
            {
                int next = j - 1;
                int nextnext = j -2;

                if (next < 0)
                {
                    next += 3;
                }

                if (nextnext < 0)
                {
                    nextnext += 3;
                }

				cornerArea_[i][j] = scale * (bcw[next] +
				                             bcw[nextnext]);
            }
		}

		triangleArea_[t[0]] += cornerArea_[i][0];
		triangleArea_[t[1]] += cornerArea_[i][1];
		triangleArea_[t[2]] += cornerArea_[i][2];
	}
}

void Mesh::update()
{
    for (int i=0;i<triangles_.size();i++)
    {
        auto& t = triangles_[i];

        for (int j=0;j<3;j++)
        {
            int index = t[j];
            t.vertices_[j] = vertices_[index];
        }
    }
}

void Mesh::CalcVertexNormals()
{
    MeshTools::CalculateTriangleAreasAndFaceNormals(*this, triangleArea_, facetNormals_);

    // calculate vertex normals
    std::vector<Real3> VertNorms(getNumVertices(), {{0,0,0}});

    // calculate the normals of each vertices
    for (int i=0;i<triangles_.size();i++){
        auto& t = triangles_[i];
        Real3 normals = facetNormals_[i];

        for (int j=0;j<3;j++){
            for (int k=0;k<3;k++){
                VertNorms[t[j]][k] += normals[k];
            }
        }
    }

    #pragma omp parallel for
    for (int i=0;i<getNumVertices();i++)
    {
        Real norm = LinAlg3x3::norm(VertNorms[i]);
        vertices_[i].normals_ = VertNorms[i]/norm;
    }
}

void Mesh::CalcVertexNormalsAreaWeighted()
{
    // calculate vertex normals
    vertexNormals_.resize(vertices_.size());
    Real3 zeroArr = {{0,0,0}};
    std::fill(vertexNormals_.begin(), vertexNormals_.end(), zeroArr);

    // calculate the normals of each vertices
    for (int i=0;i<triangles_.size();i++){
        auto& t = triangles_[i];
        Real area = triangleArea_[i];
        Real3 normals = facetNormals_[i];

        for (int j=0;j<3;j++){
            auto& vNorm = vertexNormals_[t[j]];
            vNorm = vNorm + area * normals;
        }
    }


    for (int i=0;i<vertexNormals_.size();i++){
        Real norm = LinAlg3x3::norm(vertexNormals_[i]);
        vertexNormals_[i] = vertexNormals_[i] / norm;
    }

    for (int i=0;i<vertices_.size();i++){
        vertices_[i].normals_ = vertexNormals_[i];
    }
}

std::vector<Mesh::Real3> Mesh::getVertexPositions()
{
    std::vector<Real3> v;
    for (int i=0;i<vertices_.size();i++)
    {
        v.push_back(vertices_[i].position_);
    }
    return v;
}

std::vector<Mesh::INT3> Mesh::getFaces(){
    std::vector<INT3> faces;
    for (int i=0;i<triangles_.size();i++){
        faces.push_back(triangles_[i].triangleindices_);
    }

    return faces;
}

                                    /**************************************
                                     ************   MeshTools *************
                                     *************************************/       
void MeshTools::writePLY(std::string filename, const Mesh& mesh)
{
    std::vector<Real3> verts;
    std::vector<INT3> faces;
    std::vector<Real3> normals;

    const auto& v = mesh.getvertices();
    const auto& f = mesh.gettriangles();

    verts.resize(v.size());
    faces.resize(f.size());
    normals.resize(v.size());

    for (int i=0;i<v.size();i++){
        verts[i] = v[i].position_;
        normals[i] = v[i].normals_;
    }

    for (int i=0;i<f.size();i++){
        faces[i] = f[i].triangleindices_;
    }

    writePLY(filename, verts, faces, normals);
}

void MeshTools::writePLY(std::string filename, const std::vector<Real3>& vertices, const std::vector<INT3>& faces, const std::vector<Real3>& normals)
{
    std::ofstream ofs;
    ofs.open(filename);

    ofs << "ply" << "\n";
    ofs << "format ascii 1.0\n";
    ofs << "comment Created by Yusheng Cai\n";

    int sizeVertex = vertices.size();
    int sizetriangle = faces.size();

    ofs << "element vertex " << sizeVertex << std::endl; 
    ofs << "property float x\n";
    ofs << "property float y\n";
    ofs << "property float z\n";
    ofs << "property float nx\n";
    ofs << "property float ny\n";
    ofs << "property float nz\n";
    ofs << "element face " << sizetriangle << "\n";
    ofs << "property list uchar uint vertex_indices\n";
    ofs << "end_header\n";

    ofs << std::fixed << std::setprecision(6);
    for (int i=0;i<sizeVertex;i++){
        for (int j=0;j<3;j++){
            ofs << vertices[i][j] << " ";
        }

        for (int j=0;j<3;j++){
            ofs << normals[i][j] << " ";
        }

        ofs << "\n";
    }

    for (int i=0;i<sizetriangle;i++){
        ofs << 3 << " ";
        auto& t = faces[i];

        for (int j=0;j<3;j++){
            ofs << t[j] << " ";
        }
        ofs << "\n";
    }

    ofs.close();
}


void MeshTools::writePLY(std::string filename, const std::vector<Real3>& vertices, const std::vector<INT3>& faces, \
Real factor)
{
    std::ofstream ofs;
    ofs.open(filename);

    ofs << "ply" << "\n";
    ofs << "format ascii 1.0\n";
    ofs << "comment Created by Yusheng Cai\n";

    int sizeVertex = vertices.size();
    int sizetriangle = faces.size();

    ofs << "element vertex " << sizeVertex << std::endl; 
    ofs << "property float x\n";
    ofs << "property float y\n";
    ofs << "property float z\n";
    ofs << "element face " << sizetriangle << "\n";
    ofs << "property list uchar uint vertex_indices\n";
    ofs << "end_header\n";

    ofs << std::fixed << std::setprecision(6);
    for (int i=0;i<sizeVertex;i++)
    {
        for (int j=0;j<3;j++)
        {
            ofs << vertices[i][j] * factor << " ";
        }

        ofs << "\n";
    }

    for (int i=0;i<sizetriangle;i++)
    {
        ofs << 3 << " ";
        auto& t = faces[i];

        for (int j=0;j<3;j++)
        {
            ofs << t[j] << " ";
        }
        ofs << "\n";
    }

    ofs.close();
}

void MeshTools::writePLYRGB(std::string filename, const std::vector<Real3>& vertices, const std::vector<INT3>& faces, const std::vector<Real3>& RGB)
{
    std::ofstream ofs;
    ofs.open(filename);

    ofs << "ply" << "\n";
    ofs << "format ascii 1.0\n";
    ofs << "comment Created by Yusheng Cai\n";

    int sizeVertex = vertices.size();
    int sizetriangle = faces.size();
    int sizeRGB    = RGB.size();

    ofs << "element vertex " << sizeVertex << std::endl; 
    ofs << "property float x\n";
    ofs << "property float y\n";
    ofs << "property float z\n";
    ofs << "property uchar red\n";
    ofs << "property uchar green\n";
    ofs << "property uchar blue\n";
    ofs << "element face " << sizetriangle << "\n";
    ofs << "property list uchar uint vertex_indices\n";
    ofs << "end_header\n";

    ASSERT((sizeRGB == sizeVertex), "The size of RGB value provided must agree with the size of vertex.");

    ofs << std::fixed << std::setprecision(6);
    for (int i=0;i<sizeVertex;i++)
    {
        for (int j=0;j<3;j++)
        {
            ofs << vertices[i][j] << " ";
        }

        for (int j=0;j<3;j++)
        {
            ofs << (int)RGB[i][j] << " ";
        }

        ofs << "\n";
    }

    for (int i=0;i<sizetriangle;i++)
    {
        ofs << 3 << " ";
        auto& t = faces[i];

        for (int j=0;j<3;j++)
        {
            ofs << t[j] << " ";
        }
        ofs << "\n";
    }

    ofs.close();
}

void MeshTools::writeSTL(std::string name, const std::vector<Real3>& vertices, const std::vector<INT3>& faces){

}

void MeshTools::writeSTL(std::string name, Mesh& mesh){
    std::ofstream ofs;
    ofs.open(name);

    const auto& f = mesh.gettriangles();
    const auto& n = mesh.getFaceNormals();
    
    ofs << "solid " << name << "\n";

    for (int i=0;i<f.size();i++)
    {
        ofs << "facet normal " << n[i][0] << " " << n[i][1] << " " << n[i][2] << "\n";

        ofs << "\touter loop\n";
        auto& t = f[i];

        for (int j=0;j<3;j++)
        {
            auto& pos = t.vertices_[j].position_;
            ofs << "vertex " << pos[0] << " " << pos[1] << " " << pos[2] << "\n";
        }

        ofs << "\tendloop\n";
        ofs << "endfacet\n";
    }
    ofs << "endsolid " << name;

}

bool MeshTools::readPLYlibr(std::string& filename, Mesh& mesh)
{
    happly::PLYData plydata(filename);
    std::vector<std::array<double,3>> vPos = plydata.getVertexPositions();
    std::vector<std::vector<size_t>> fInd = plydata.getFaceIndices<size_t>();
    auto names = plydata.getElement("vertex").getPropertyNames();

    // first we update the mesh
    auto& vertices = mesh.accessvertices();
    auto& triangles= mesh.accesstriangles();
    auto& verticesPos = mesh.accessverticesPos();

    vertices.clear();
    triangles.clear();

    vertices.resize(vPos.size());
    verticesPos.resize(vPos.size());
    triangles.resize(fInd.size());

    for (int i=0;i<vertices.size();i++){
        auto& v = vertices[i];
        for (int j=0;j<3;j++){
            v.position_[j] = vPos[i][j];
        }
    }

    for (int i=0;i<verticesPos.size();i++){
        auto& y = verticesPos[i];
        for (int j=0;j<3;j++){
            y[j] = vPos[i][j];
        }
    }


    for (int i=0;i<triangles.size();i++){
        auto& t = triangles[i];

        ASSERT((fInd[i].size() == 3), "We are reading triangles here while the provided face is " << fInd[i].size());

        for (int j=0;j<3;j++){
            t[j] = fInd[i][j];
            t.vertices_[j] = vertices[fInd[i][j]];
        }
    }

    mesh.CalcVertexNormals();
    mesh.update();

    return true;
}

bool MeshTools::readPLY(std::string& filename, Mesh& mesh)
{
    auto& vertices = mesh.accessvertices();
    auto& triangles= mesh.accesstriangles();
    vertices.clear();
    triangles.clear();

    // open the file
    std::ifstream ifs;
    std::stringstream ss;
    ifs.open(filename);

    if (! ifs.is_open()){
        return false;
    }

    std::string sentence;
    int numfaces;
    int numvertex;
    while (std::getline(ifs, sentence)){
        ss.str(sentence);
        std::string token;

        std::vector<std::string> vectorstring;
        while (ss >> token){
            vectorstring.push_back(token);
        }

        if (vectorstring[0] == "end_header"){
            break;
        }

        if (vectorstring[0] == "format"){
            ASSERT((vectorstring[1] == "ascii"), "Currently do not support non-ascii format ply files.");
        }

        if (vectorstring[0] == "element"){
            ASSERT((vectorstring.size() == 3), "The element can only have 3.");
            if (vectorstring[1] == "vertex"){
                numvertex = StringTools::StringToType<int>(vectorstring[2]);
            }

            if (vectorstring[1] == "face"){
                numfaces = StringTools::StringToType<int>(vectorstring[2]);
            }
        }

        ss.clear();
    }

    ss.clear();
    // read in the vertices as well as their normals     
    std::string datasentence;
    for (int i=0;i<numvertex;i++)
    {
        std::getline(ifs, datasentence);
        ss.str(datasentence);
        std::string data;
        std::vector<std::string> vectordata;

        vertex v;

        while (ss >> data)
        {
            vectordata.push_back(data);
        }

        ASSERT((vectordata.size() == 6), "The ply file must contain x, y, z, nx, ny ,nz");

        Real3 position;
        Real3 normals;

        for (int j=0;j<3;j++)
        {
            position[j] = StringTools::StringToType<Real>(vectordata[j]);
            normals[j] = StringTools::StringToType<Real>(vectordata[j+3]);
        }

        v.position_ = position;
        v.normals_ = normals;

        vertices.push_back(v);

        ss.clear();
    }

    ss.clear();
    // read in the triangles
    std::string trianglesentence_;
    for (int i=0;i<numfaces;i++)
    {
        std::getline(ifs,trianglesentence_);
        ss.str(trianglesentence_);
        std::string data;
        std::vector<std::string> vectordata;

        while (ss >> data){
            vectordata.push_back(data);
        }

        ASSERT((vectordata.size() == 4), "The triangle file must contain 3 f1 f2 f3");

        triangle t;

        INT3 faceid;

        for (int j=0;j<3;j++){faceid[j] = StringTools::StringToType<int>(vectordata[j+1]);}

        t.triangleindices_ = faceid;
        for (int j=0;j<3;j++){
            t.vertices_[j] = vertices[faceid[j]];
        }

        triangles.push_back(t);
        ss.clear();
    }

    ifs.close();

    return true;
}

MeshTools::Real3 MeshTools::calculateShift(const Real3& vec1, const Real3& vec2, const Real3& boxLength)
{
    Real3 diff;
    for (int i=0;i<3;i++)
    {
        Real d = vec1[i] - vec2[i];

        if (d > (0.5 * boxLength[i])) diff[i] = -boxLength[i];
        else if (d < (- 0.5 * boxLength[i])) diff[i] = boxLength[i];
        else diff[i] = 0.0;
    }
    
    return diff;
}

void MeshTools::calculateDistance(const Real3& vec1, const Real3& vec2, const Real3& boxLength, Real3& distance, Real& distsq)
{
    distsq = 0.0;
    for (int i=0;i<3;i++)
    {
        Real d = vec1[i] - vec2[i];

        if (d > (0.5 * boxLength[i])) distance[i] = d-boxLength[i];
        else if (d < (- 0.5 * boxLength[i])) distance[i] = d+boxLength[i];
        else distance[i] = d;

        distsq += distance[i] * distance[i];
    }
}

bool MeshTools::isPeriodicEdge(const Real3& vec1, const Real3& vec2, Real3& newarr, const Real3& boxLength)
{
    // calculate pbc shift b/t vec1 and vec2 in some box with boxlength
    Real3 diff = calculateShift(vec1, vec2, boxLength);

    // shifted vector 
    newarr = vec1 + diff;

    // check if this is periodic edge
    for (int i=0;i<3;i++){
        if (std::abs(diff[i]) > (boxLength[i] * 0.5)){
            return true;
        }
    }

    return false;
}

void MeshTools::MapVerticesToFaces(const Mesh& mesh, std::vector<std::vector<int>>& map){
    const auto& vertices = mesh.getvertices();
    const auto& triangles= mesh.gettriangles();

    map.clear();
    map.resize(vertices.size());

    for (int i=0;i<triangles.size();i++){
        auto t = triangles[i].triangleindices_;
        for (int j=0;j<3;j++){
            map[t[j]].push_back(i);
        }
    }
}

void MeshTools::CalculateTriangleAreasAndFaceNormals(Mesh& mesh, std::vector<Real>& Areas, std::vector<Real3>& Normals)
{
    auto& triangles = mesh.accesstriangles();
    auto& vertices  = mesh.accessvertices();

    Areas.clear();
    Areas.resize(triangles.size());

    Normals.clear();
    Normals.resize(triangles.size());

    // calculate the area of the triangles as well as the normals of the faces of triangles
    #pragma omp parallel for
    for (int i=0;i<triangles.size();i++){
        auto& t = triangles[i];

        int index1 = t[0];
        int index2 = t[1];
        int index3 = t[2];

        Real3 diff1 = {};
        Real3 diff2 = {};
        Real norm1, norm2;

        mesh.getVertexDistance(vertices[index2], vertices[index1], diff1, norm1);
        mesh.getVertexDistance(vertices[index3], vertices[index1], diff2, norm2);

        Real3 crossProduct = LinAlg3x3::CrossProduct(diff1, diff2);
        Real norm = LinAlg3x3::norm(crossProduct);
        Real a    = norm*0.5;
        Real3 n   = crossProduct/norm;

        Normals[i] = n;
        Areas[i]   = a;
    }
}

void MeshTools::CalculateVertexNeighbors(const Mesh& mesh, std::vector<std::vector<int>>& neighborIndices)
{
    auto& vertices = mesh.getvertices();
    auto& triangles= mesh.gettriangles();

    neighborIndices.clear();
    neighborIndices.resize(vertices.size());

    for (int i=0;i<triangles.size();i++){
        auto& t = triangles[i];

        for (int j=0;j<3;j++){
            int index1 = t[j];
            for (int k=0;k<3;k++){
                if (j != k){
                    int index2 = t[k];
                    auto it = std::find(neighborIndices[index1].begin(), neighborIndices[index1].end(), index2);
                    if (it == neighborIndices[index1].end()){
                        neighborIndices[index1].push_back(index2);
                    }
                }
            }
        }
    }
}

void MeshTools::MapEdgeToFace(const Mesh& mesh, std::map<INT2, std::vector<int>>& mapEdgeToFace, bool assert){
    std::vector<std::vector<INT2>> temp;
    MapEdgeToFace(mesh, mapEdgeToFace, temp, assert);
}

void MeshTools::MapEdgeToFace(const Mesh& mesh, std::map<INT2, std::vector<int>>& mapEdgeToFace, std::vector<std::vector<INT2>>& mapVertexToEdge, bool assert)
{
    // get the triangles and vertices from mesh object
    const auto& triangles = mesh.gettriangles();
    const auto& vertices  = mesh.getvertices();

    // clear the inputs 
    mapVertexToEdge.clear();
    mapVertexToEdge.resize(vertices.size());
    mapEdgeToFace.clear();

    // iterate over all the triangles
    for (int i=0;i<triangles.size();i++){
        // find the triangles indices
        auto& t = triangles[i];

        for (int j=0;j<3;j++){
            // get the 2 adjacent indices
            int idx1 = t[j];
            int idx2 = t[(j+1)%3];

            // find the min/max of the 2
            int minIndex = std::min(idx1, idx2);
            int maxIndex = std::max(idx1, idx2);
            INT2 arr = {{minIndex, maxIndex}};

            // insert into map
            auto it = mapEdgeToFace.find(arr);
            if (it == mapEdgeToFace.end()){
                std::vector<int> faceIndices;
                faceIndices.push_back(i);

                mapEdgeToFace.insert(std::make_pair(arr,faceIndices));
            }
            else{
                it -> second.push_back(i);
            }
        }

        // map vertex to edges 
        for (int j=0;j<3;j++){
            int idx1 = t[j];
            int idx2 = t[(j+1)%3];
            int minIdx = std::min(idx1, idx2);
            int maxIdx = std::max(idx1, idx2);

            INT2 indices = {{minIdx, maxIdx}};

            for (int k=0;k<2;k++){
                // find if the edge is already in the vector
                auto f = std::find(mapVertexToEdge[indices[k]].begin(), mapVertexToEdge[indices[k]].end(), indices);

                if (f == mapVertexToEdge[indices[k]].end()){
                    mapVertexToEdge[indices[k]].push_back(indices);
                }
            }
        }
    }

    // do a check to make sure that each edge is shared by at least 1 and at most 2 faces
    if (assert){
        int id = 0;
        for (auto it=mapEdgeToFace.begin(); it != mapEdgeToFace.end();it++){
            if (it -> second.size() > 2){
                std::cout << "This is for edge " << id << " consisting of vertex " << it ->first <<"\n";
                for (auto s : it -> second){
                    std::cout << "Face it corresponds to : " << s <<"\n";
                }
            }
            ASSERT((it->second.size() == 1 || it -> second.size() ==2), "The number of faces shared by an edge needs to be either 1 or 2 \
            , however the calculated result shows " << it -> second.size());
            id ++;
        }
    }

    return;
}

void MeshTools::CalculateBoundaryVertices(const Mesh& mesh, const std::map<INT2, std::vector<int>>& mapEdgeToFace, std::vector<bool>& boundaryIndicator)
{
    const auto& vertices = mesh.getvertices();

    boundaryIndicator.clear();
    boundaryIndicator.resize(vertices.size(), false);

    for (auto it = mapEdgeToFace.begin(); it != mapEdgeToFace.end(); it ++){
        int size = it -> second.size();
        INT2 arr = it -> first;

        // if the edge size == 1, then it is a boundary edge
        if (size == 1){
            boundaryIndicator[arr[0]] = true;
            boundaryIndicator[arr[1]] = true;
        }
    }
}

void MeshTools::CalculateBoundaryVertices(const Mesh& mesh, std::vector<bool>& boundaryIndicator)
{
    std::map<INT2, std::vector<int>> mapEdgeToFace;

    MapEdgeToFace(mesh, mapEdgeToFace);

    CalculateBoundaryVertices(mesh, mapEdgeToFace, boundaryIndicator);
}

bool MeshTools::IsBoundary(int index, const std::vector<bool>& boundaryIndicator)
{
    return boundaryIndicator[index];
}

void MeshTools::ConvertToNonPBCMesh(Mesh& mesh, bool AddNewTriangles){
    std::vector<Real3> newv;
    std::vector<INT3> newf;
    ConvertToNonPBCMesh(mesh, newv, newf, AddNewTriangles);

    auto& oldv = mesh.accessvertices();
    auto& oldf = mesh.accesstriangles();
    oldv.clear(); oldv.resize(newv.size());
    oldf.clear(); oldf.resize(newf.size());

    for (int i=0;i<newv.size();i++){
        oldv[i].position_ = newv[i];
    }

    for (int i=0;i<newf.size();i++){
        oldf[i].triangleindices_ = newf[i];
    }
}

void MeshTools::ConvertToNonPBCMesh(Mesh& mesh, std::vector<Real3>& vertices, std::vector<INT3>& triangles, bool AddNewTriangles)
{
    // if periodic, then do something, else do nothing 
    if (mesh.isPeriodic()){
        auto MeshVertices = mesh.getvertices();
        auto MeshTriangles = mesh.gettriangles();

        // first let's copy all the vertices 
        vertices.clear();
        for (auto& v : MeshVertices){
            vertices.push_back(v.position_);
        }

        std::map<int,int> MapShiftedIndicesToOld;

        // check if triangle is periodic 
        for (auto& t : MeshTriangles){
            // find all the edge lenghs
            std::map<INT2, bool> mapEdge;
            bool periodicTriangle = MeshTools::IsPeriodicTriangle(MeshVertices, t.triangleindices_, mesh.getBoxLength(), mapEdge);

            // if this particular triangle is not periodic 
            if (! periodicTriangle){
                triangles.push_back(t.triangleindices_);
            }
            // if it's periodic triangle, then we push back 3 new vertices 
            else{
                if (AddNewTriangles){
                    int shift_index, newIndex;
                    Real3 verticesNew, verticesDiff;
                    Real distsq;

                    if (!MeshTools::FindPeriodicShiftIndex(mapEdge, shift_index)){shift_index=0;std::cout << "Shift index is not found." << std::endl;}

                    if (!Algorithm::FindInMap(MapShiftedIndicesToOld, t[shift_index], newIndex)){
                        int idx1=t[shift_index];
                        int idx2=t[(shift_index+1)%3];

                        mesh.getVertexDistance(MeshVertices[idx1].position_, MeshVertices[idx2].position_, verticesDiff, distsq);
                        verticesNew = MeshVertices[idx2].position_ + verticesDiff;

                        newIndex = vertices.size();
                        vertices.push_back(verticesNew);

                        Algorithm::InsertInMap(MapShiftedIndicesToOld, t[shift_index], newIndex);
                    }

                    INT3 newt;
                    newt[shift_index] = newIndex;
                    newt[(shift_index+1)%3] = t[(shift_index + 1)%3]; 
                    newt[(shift_index+2)%3] = t[(shift_index + 2)%3];

                    triangles.push_back(newt);
                }
            }
        }
    }
}

bool MeshTools::IsPeriodicTriangle(std::vector<vertex>& Vertices,INT3& face, Real3 BoxLength)
{
    bool IsPeriodic=false;
    for (int i=0;i<3;i++){
        int index1 = face[i];
        int index2 = face[(i+1) % 3]; 
        Real3 diff;
        for (int j=0;j<3;j++){
            diff[j] = Vertices[index1].position_[j] - Vertices[index2].position_[j];

            if (BoxLength[j] != 0){
                if (std::abs(diff[j]) >= 0.5 * BoxLength[j]){
                    return true;
                }
            }
        }
    }

    return false;
}

bool MeshTools::IsPeriodicTriangle(std::vector<vertex>& Vertices,INT3& face, Real3 BoxLength, std::map<INT2,bool>& mapEdge)
{
    mapEdge.clear();
    for (int i=0;i<3;i++){
        int index1 = face[i];
        int index2 = face[(i+1) % 3]; 
        Real3 diff;
        bool isPeriodicEdge=false;
        for (int j=0;j<3;j++){
            diff[j] = Vertices[index1].position_[j] - Vertices[index2].position_[j];

            if (BoxLength[j] != 0){
                if (std::abs(diff[j]) >= 0.5 * BoxLength[j]){
                    isPeriodicEdge=true;
                    break;
                }
            }
        }
        INT2 edge = MeshTools::makeEdge(i, (i+1)%3);
        Algorithm::InsertInMap(mapEdge, edge, isPeriodicEdge);
    }

    for (auto it=mapEdge.begin(); it != mapEdge.end(); it++){
        if (it ->second == true){
            return true;
        }
    }

    return false;
}



bool MeshTools::IsPeriodicTriangle(const Mesh& mesh, int faceIndex, std::map<INT2,bool>& mapEdge){
    const auto& f = mesh.gettriangles();
    const auto& v = mesh.getvertices();
    Real3 boxLength = mesh.getBoxLength();
    mapEdge.clear();

    INT3 t = f[faceIndex].triangleindices_;
    for (int i=0;i<3;i++){
        int index1 = t[i];
        int index2 = t[(i+1) % 3];
        Real3 diff;
        for (int j=0;j<3;j++){
            diff[j] = v[index1].position_[j] - v[index2].position_[j];

            if (boxLength[j] != 0){
                if (std::abs(diff[j]) > 0.5 * boxLength[j]){
                    INT2 edge = MeshTools::makeEdge(i, (i+1) % 3);
                    Algorithm::InsertInMap(mapEdge, edge, true);
                    break;
                }
            }
        }
    }

    for (auto it=mapEdge.begin(); it != mapEdge.end(); it++){
        if (it ->second == true){
            return true;
        }
    }

    return false;
}

bool MeshTools::IsPeriodicTriangle(const Mesh& mesh, int faceIndex){
    const auto& f = mesh.gettriangles();
    const auto& v = mesh.getvertices();
    Real3 boxLength = mesh.getBoxLength();
    INT3 t = f[faceIndex].triangleindices_;
    for (int i=0;i<3;i++){
        int index1 = t[i];
        int index2 = t[(i+1) % 3];
        Real3 diff;
        for (int j=0;j<3;j++){
            diff[j] = v[index1].position_[j] - v[index2].position_[j];

            if (boxLength[j] != 0){
                if (std::abs(diff[j] >= 0.5 * boxLength[j])){return true;}
            }
        }
    }

    return false;
}

bool MeshTools::FindPeriodicShiftIndex(std::map<INT2,bool>& mapEdge, int& index){
    int ind=0;
    INT2 nonPeriodicEdge;
    bool hasNonPeriodicEdge=false;
    for (auto it = mapEdge.begin(); it != mapEdge.end(); it ++){
        // if edge is not periodic
        if (it ->second == false){
            hasNonPeriodicEdge=true;
            nonPeriodicEdge = it->first;
            break;
        }
    }

    if (! hasNonPeriodicEdge){
        return false;
    }
    else{
        for (int i=0;i<3;i++){
            if (i!=nonPeriodicEdge[0] && i != nonPeriodicEdge[1]){
                index = i;
                return true;
            }
        }
    }
}

void MeshTools::ShiftPeriodicTriangle(const std::vector<vertex>& Vertices, const INT3& face, Real3 BoxLength, Real3& A, Real3& B, Real3& C)
{
    // half box
    Real3 half_box = BoxLength * 0.5;

    // first shift B C wrt to A
    A = Vertices[face[0]].position_;
    B = Vertices[face[1]].position_;
    C = Vertices[face[2]].position_;
    Real3 shiftB = MeshTools::calculateShift(B, A, BoxLength);
    Real3 shiftC = MeshTools::calculateShift(C, A, BoxLength);
    B = B + shiftB;
    C = C + shiftC;

    // calculate the triangle center
    Real3 center = (A + B + C) * 1.0/3.0;

    // calculate shift
    Real3 shiftTriangle = MeshTools::calculateShift(center, half_box, BoxLength);

    A = A + shiftTriangle;
    B = B + shiftTriangle;
    C = C + shiftTriangle;
}

void MeshTools::ShiftPeriodicTriangle(const std::vector<vertex>& Vertices, const INT3& face, Real3 BoxLength, const Real3& point, Real3& A, Real3& B, Real3& C){
    // half box
    Real3 half_box = BoxLength * 0.5;

    // first shift B C wrt to A
    A = Vertices[face[0]].position_;
    B = Vertices[face[1]].position_;
    C = Vertices[face[2]].position_;
    Real3 shiftA = MeshTools::calculateShift(A, point, BoxLength);
    Real3 shiftB = MeshTools::calculateShift(B, point, BoxLength);
    Real3 shiftC = MeshTools::calculateShift(C, point, BoxLength);

    A = A + shiftA; B = B + shiftB; C = C+ shiftC;
}

MeshTools::INT2 MeshTools::makeEdge(int i, int j)
{
    int minIndex = std::min(i,j);
    int maxIndex = std::max(i,j);

    INT2 ret = {{minIndex, maxIndex}};

    return ret;
}

void MeshTools::MapEdgeToOpposingVertices(Mesh& mesh, std::map<INT2, std::vector<int>>& mapEdgeToFace, std::map<INT2, std::vector<int>>& MapEdgeToOppoVertices)
{
    const auto& tri = mesh.gettriangles();

    // clear the map at which we want to write to  
    MapEdgeToOppoVertices.clear();

    // iterate over the edges 
    for (auto it = mapEdgeToFace.begin(); it != mapEdgeToFace.end(); it ++){
        std::vector<int> faces = it -> second;
        INT2 edge = it -> first;
        std::vector<int> opposingPoints;

        int numf = faces.size();

        // only perform calculation if we are working with non boundary edge --> (edges shared by 2 faces)
        // iterate over the 2 faces 
        for (int f : faces){
            auto& TriIndices = tri[f].triangleindices_;

            // iterate over the triangular indices of each of the face 
            for (int id : TriIndices){
                bool IsInEdge = std::find(edge.begin(), edge.end(), id) != edge.end();

                if ( ! IsInEdge){
                    opposingPoints.push_back(id);
                }
            }
        }

        MapEdgeToOppoVertices.insert(std::make_pair(edge, opposingPoints));
    }
}

void MeshTools::CutMesh(Mesh& mesh, std::vector<INT3>& faces, std::vector<Real3>& vertices, Real3 volume, bool below)
{
    vertices.clear();
    faces.clear();

    const auto& verts = mesh.getvertices();
    const auto& tri   = mesh.gettriangles();

    int index=0;
    std::map<int,int> MapOldIndexToNew;
    for (int i=0;i<verts.size();i++){
        auto& v = verts[i];
        if (! below){
            if (v.position_[0] >= volume[0] && v.position_[1] >= volume[1] && v.position_[2] >= volume[2]){
                int newindex = index;
                vertices.push_back(v.position_);
                MapOldIndexToNew.insert(std::make_pair(i, newindex));
                index ++;
            }
        }
        else{
            if (v.position_[0] <= volume[0] && v.position_[1] <= volume[1] && v.position_[2] <= volume[2]){
                int newindex = index;
                vertices.push_back(v.position_);
                MapOldIndexToNew.insert(std::make_pair(i, newindex));
                index ++;
            }
        }
    }

    for (auto& t : tri){
        bool it1 = MapOldIndexToNew.find(t.triangleindices_[0]) != MapOldIndexToNew.end();
        bool it2 = MapOldIndexToNew.find(t.triangleindices_[1]) != MapOldIndexToNew.end();
        bool it3 = MapOldIndexToNew.find(t.triangleindices_[2]) != MapOldIndexToNew.end();

        if (it1 && it2 && it3){
            INT3 newt;
            for (int i=0;i<3;i++){
                auto it = MapOldIndexToNew.find(t.triangleindices_[i]);
                int newIndex = it -> second;
                newt[i] = newIndex;
            }
            faces.push_back(newt);
        }
    }
}

void MeshTools::CutMesh(Mesh& mesh, Real3& volume, bool below){
    const auto& verts = mesh.getvertices();
    const auto& tri   = mesh.gettriangles();
    std::vector<vertex> newV;
    std::vector<triangle> newT;

    int index=0;
    std::map<int,int> MapOldIndexToNew;
    for (int i=0;i<verts.size();i++){
        auto& v = verts[i];
        if (! below){
            if (v.position_[0] >= volume[0] && v.position_[1] >= volume[1] && v.position_[2] >= volume[2]){
                newV.push_back(v);
                MapOldIndexToNew.insert(std::make_pair(i, index));
                index ++;
            }
        }
        else{
            if (v.position_[0] <= volume[0] && v.position_[1] <= volume[1] && v.position_[2] <= volume[2]){
                newV.push_back(v);
                MapOldIndexToNew.insert(std::make_pair(i, index));
                index ++;
            }
        }
    }

    for (auto& t : tri){
        bool it1 = MapOldIndexToNew.find(t.triangleindices_[0]) != MapOldIndexToNew.end();
        bool it2 = MapOldIndexToNew.find(t.triangleindices_[1]) != MapOldIndexToNew.end();
        bool it3 = MapOldIndexToNew.find(t.triangleindices_[2]) != MapOldIndexToNew.end();

        if (it1 && it2 && it3){
            triangle newt;
            for (int i=0;i<3;i++){
                newt.triangleindices_[0] = MapOldIndexToNew[t.triangleindices_[0]];
                newt.triangleindices_[1] = MapOldIndexToNew[t.triangleindices_[1]];
                newt.triangleindices_[2] = MapOldIndexToNew[t.triangleindices_[2]];
            }
            newT.push_back(newt);
        }
    }

    mesh.SetVerticesAndTriangles(newV, newT);
}

void MeshTools::CalculateCornerArea(Mesh& mesh, std::vector<Real3>& CornerArea, std::vector<Real>& VertexArea)
{
    const auto& vertices = mesh.getvertices();
    const auto& triangles= mesh.gettriangles();
    int nf = triangles.size();
    int nv = vertices.size();

    // resize corner area
    CornerArea.clear();
    CornerArea.resize(nf);

    // resize vertex area
    VertexArea.clear();
    VertexArea.resize(nv,0.0);

	for (int i = 0; i < nf; i++) {
		// Edges
        auto& t = triangles[i];
        int index1 = t[0];
        int index2 = t[1];
        int index3 = t[2];

        Real3 edge1, edge2, edge3;
        Real edge1sq, edge2sq, edge3sq;

        mesh.getVertexDistance(vertices[index3], vertices[index2], edge1, edge1sq);
        mesh.getVertexDistance(vertices[index1], vertices[index3], edge2, edge2sq);
        mesh.getVertexDistance(vertices[index2], vertices[index1], edge3, edge3sq);

		// Compute corner weights
        Real3 cross = LinAlg3x3::CrossProduct(edge1, edge2);
        Real  area  = 0.5 * LinAlg3x3::norm(cross);

        Real e1 = LinAlg3x3::norm(edge1);
        Real e2 = LinAlg3x3::norm(edge2);
        Real e3 = LinAlg3x3::norm(edge3);

		Real3 l2 = { e1*e1, e2*e2, e3*e3};

		// Barycentric weights of circumcenter
		Real3 bcw = { l2[0] * (l2[1] + l2[2] - l2[0]),
		                 l2[1] * (l2[2] + l2[0] - l2[1]),
		                 l2[2] * (l2[0] + l2[1] - l2[2]) };

		if (bcw[0] <= 0.0f) {
			CornerArea[i][1] = -0.25f * l2[2] * area/LinAlg3x3::DotProduct(edge1, edge3);
			CornerArea[i][2] = -0.25f * l2[1] * area/LinAlg3x3::DotProduct(edge1, edge2);
			CornerArea[i][0] = area - CornerArea[i][1] - CornerArea[i][2];
		} else if (bcw[1] <= 0.0f) {
			CornerArea[i][2] = -0.25f * l2[0] * area/LinAlg3x3::DotProduct(edge2, edge1);
			CornerArea[i][0] = -0.25f * l2[2] * area/LinAlg3x3::DotProduct(edge2, edge3);
			CornerArea[i][1] = area - CornerArea[i][2] - CornerArea[i][0];
		} else if (bcw[2] <= 0.0f) {
			CornerArea[i][0] = -0.25f * l2[1] * area/LinAlg3x3::DotProduct(edge3, edge2);
			CornerArea[i][1] = -0.25f * l2[0] * area/LinAlg3x3::DotProduct(edge3, edge1);
			CornerArea[i][2] = area - CornerArea[i][0] - CornerArea[i][1];
		} else {
			float scale = 0.5f * area / (bcw[0] + bcw[1] + bcw[2]);
			for (int j = 0; j < 3; j++){
                int next = j - 1;
                int nextnext = j -2;

                if (next < 0){
                    next += 3;
                }

                if (nextnext < 0){
                    nextnext += 3;
                }

				CornerArea[i][j] = scale * (bcw[next] +
				                             bcw[nextnext]);
            }
		}

		VertexArea[t[0]] += CornerArea[i][0];
		VertexArea[t[1]] += CornerArea[i][1];
		VertexArea[t[2]] += CornerArea[i][2];
	}
}

bool MeshTools::MTRayTriangleIntersection(Real3& A, Real3& B, Real3& C, Real3& O, Real3& D, Real& t, Real& u, Real& v)
{
    Real epsilon = 1e-8;
    
    // declare the variables 
    Real3 T, E1, E2, P, Q;
    Real det, invdet;

    // find E1=B-A and E2=C-A
    E1 = B - A;
    E2 = C - A;

    P = LinAlg3x3::CrossProduct(D, E2);
    det= LinAlg3x3::DotProduct(P, E1);

    // if det is close to 0, then the ray and triangle are parallel
    if (std::abs(det) < epsilon){
        return false;
    }

    // calculate inverse of det  
    invdet = 1.0 / det;
    T = O - A;

    // find u
    u = LinAlg3x3::DotProduct(T, P) * invdet;
    if (u < 0 || u > 1) return false;

    // find v
    Q = LinAlg3x3::CrossProduct(T,E1);
    v = LinAlg3x3::DotProduct(D, Q) * invdet;
    if (v < 0 || u+v > 1) return false;

    // find t
    t = LinAlg3x3::DotProduct(E2, Q) * invdet;

    return true;
}

void MeshTools::writeNonPBCMesh(std::string name, Mesh& mesh)
{
    if (mesh.isPeriodic())
    {
        std::vector<Real3> tempVertices;
        std::vector<INT3> tempFaces;
        MeshTools::ConvertToNonPBCMesh(mesh, tempVertices, tempFaces);

        MeshTools::writePLY(name, tempVertices, tempFaces);
    }
    else
    {
        std::cout << "WARNING: You asked to print NON PBC Mesh while the Mesh is not PBC." << "\n";
    }
}

void MeshTools::writeNonPeriodicTriangleIndices(std::string name, Mesh& mesh)
{
    std::ofstream ofs;
    ofs.open(name);
    std::vector<int> NonPeriodicTriangleIndices;

    const auto& f = mesh.gettriangles();
    const auto& v = mesh.getvertices();

    if (mesh.isPeriodic()){
        for (int i =0;i<f.size();i++){
            auto t = f[i].triangleindices_;
            if (! MeshTools::IsPeriodicTriangle(const_cast<std::vector<vertex>&>(v), t, mesh.getBoxLength())){
                NonPeriodicTriangleIndices.push_back(i);
            }
        }
    }
    else{
        NonPeriodicTriangleIndices.resize(f.size());
        std::iota(NonPeriodicTriangleIndices.begin(), NonPeriodicTriangleIndices.end(),0);
    }

    for (int index : NonPeriodicTriangleIndices){
        ofs << index << " ";
    }

    ofs.close();
}

void MeshTools::writeMeshArea(std::string filename, Mesh& mesh)
{
    // calculate the area and facet normals 
    std::vector<Real> triangleA;
    std::vector<Real3> facetN;
    MeshTools::CalculateTriangleAreasAndFaceNormals(mesh, triangleA, facetN);

    std::ofstream ofs;
    ofs.open(filename);

    for (int i=0;i<triangleA.size();i++){
        ofs << triangleA[i] << "\n";
    }

    ofs.close();
}

void MeshTools::writeCuttedMesh(std::string filename, Mesh& mesh, Real3& volume)
{
    std::vector<INT3> faces;
    std::vector<Real3> verts;
    CutMesh(mesh, faces, verts, volume);
    writePLY(filename, verts, faces);
}



void MeshTools::CheckDegenerateTriangle(Mesh& mesh, \
                                                    std::vector<int>& MergeFaces, \
                                                    std::vector<INT2>& MergeVertices)
{
    const auto& f = mesh.gettriangles();
    const auto& v = mesh.getvertices();
    MergeFaces.clear();
    Real epsilon=1e-8;
    Real merge_epsilon=1e-4;

    // what we check is if one side is whether or not a + b = c 
    for (int i=0;i<f.size();i++){
        INT3 indices = f[i].triangleindices_;

        // get all three vertices of the triangle
        Real3 A = v[indices[0]].position_;
        Real3 B = v[indices[1]].position_;
        Real3 C = v[indices[2]].position_;

        // obtain the side length
        Real3 vec;
        Real3 triangle_length;
        Real AB, BC, CA;
        mesh.getVertexDistance(A, B, vec, AB);
        mesh.getVertexDistance(B, C, vec, BC);
        mesh.getVertexDistance(C, A, vec, CA);

        // AB + BC > CA
        // BC + CA > AB
        // CA + AB > BC
        Real diff1 = AB + BC - CA;
        Real diff2 = BC + CA - AB;
        Real diff3 = CA + AB - BC;
 
        //ASSERT((diff1 >= 0 && diff2 >= 0 && diff3 >= 0), "Something weird with triangles.");

        if ((AB < merge_epsilon) || (BC < merge_epsilon) || (CA < merge_epsilon)){
            MergeFaces.push_back(i);
            INT2 verts_ind1, verts_ind2,  verts_ind3;

            verts_ind1 = {{indices[0], indices[1]}};
            verts_ind2 = {{indices[1], indices[2]}};
            verts_ind3 = {{indices[2], indices[0]}};
            Algorithm::sort(verts_ind1); Algorithm::sort(verts_ind2); Algorithm::sort(verts_ind3);

            if ((! Algorithm::contain(MergeVertices, verts_ind1)) && (AB<merge_epsilon)){MergeVertices.push_back(verts_ind1);}
            if ((! Algorithm::contain(MergeVertices, verts_ind2)) && (BC<merge_epsilon)){MergeVertices.push_back(verts_ind2);}
            if ((! Algorithm::contain(MergeVertices, verts_ind3)) && (CA<merge_epsilon)){MergeVertices.push_back(verts_ind3);}
        }
    }
}

bool MeshTools::decimateDegenerateTriangle(Mesh& mesh)
{
    std::vector<int> MergeFaceIndices;
    std::vector<INT2> MergeVerts;
    CheckDegenerateTriangle(mesh, MergeFaceIndices, MergeVerts);

    if (MergeFaceIndices.size() != 0){
        std::cout << "MergeFace Indices = " << MergeFaceIndices << "\n";
        auto& triangles = mesh.accesstriangles();
        auto& verts     = mesh.accessvertices();

        // declare new triangles and vertices
        std::vector<triangle> newT;
        std::vector<vertex>   newV;

        int nv = verts.size();
        int nf = triangles.size();

        std::vector<int> MapOldIndicesToNew = Algorithm::arange(0,nv,1);
        std::vector<int> Min(MergeVerts.size());
        std::vector<int> Max(MergeVerts.size());

        for (int i=0;i<MergeVerts.size();i++){
            Min[i] = MergeVerts[i][0];
            Max[i] = MergeVerts[i][1];
        }
        std::cout << "Min = " << Min << "\n";
        std::cout << "Max = " << Max << '\n';

        int index=0;
        for (int i=0;i<nv ;i++){
            int ind;
            if (Algorithm::contain(Max, i, ind)){
                MapOldIndicesToNew[i] = MapOldIndicesToNew[Min[ind]];
            }
            else{
                MapOldIndicesToNew[i] = index;
                newV.push_back(verts[i]);
                index ++;
            }
        }

        // obtain the new triangle indices 
        for (int i=0;i<triangles.size();i++){
            if (! Algorithm::contain(MergeFaceIndices,i)){
                triangle t;
                auto ind  = triangles[i].triangleindices_;
                t.triangleindices_ = {{MapOldIndicesToNew[ind[0]], MapOldIndicesToNew[ind[1]], MapOldIndicesToNew[ind[2]]}};
                if (Algorithm::is_unique(t.triangleindices_))
                {
                    newT.push_back(t);
                }
            }
        }

        triangles.clear();
        triangles = newT;
        verts.clear();
        verts = newV;

        return true;
    }

    return false;
}

void MeshTools::CorrectMesh(Mesh& mesh, std::vector<int>& FaceIndices)
{
    auto& triangles = mesh.accesstriangles();
    auto& vertices  = mesh.accessvertices();

    // declare the new triangles 
    std::vector<triangle> newT;
}

void MeshTools::CalculateBoundaryBarycenter(Mesh& mesh, std::vector<bool>& boundary_indicator, Real3& barycenter)
{
    // shift all boundary vertices wrt to the first vertex
    const auto& verts = mesh.getvertices();
    barycenter=verts[0].position_;

    for (int i=1;i<verts.size();i++){
        Real3 newpos = mesh.getShiftedVertexPosition(verts[i], verts[0]);
        barycenter = barycenter + newpos;
    }

    // average
    barycenter = barycenter * 1.0/verts.size();

    // shift into box
    Real3 shift = mesh.getShiftIntoBox(barycenter);

    barycenter = barycenter + shift;
}

bool MeshTools::IsIsolatedFace(Mesh& mesh, int faceIndex, const std::map<INT2, std::vector<int>>& mapEdgeToFace){
    const auto& faces = mesh.gettriangles();
    INT3 vInd  = faces[faceIndex].triangleindices_;

    int numBoundary=0;
    for (int i=0;i<3;i++){
        int next = (i+1) % 3;
        INT2 edge = MeshTools::makeEdge(vInd[i],vInd[next]);
        std::vector<int> faces;
        bool found = Algorithm::FindInMap(mapEdgeToFace, edge, faces);
        ASSERT((found), "The edge " << edge << " is not found.");

        if (faces.size() == 1){
            numBoundary += 1;
        }
    }

    if (numBoundary == 3){
        return true;
    }

    return false;
}

bool MeshTools::IsTeethlikeFace(Mesh& mesh, int faceIndex, const std::map<INT2, std::vector<int>>& mapEdgeToFace){
    const auto& faces = mesh.gettriangles();
    INT3 vInd  = faces[faceIndex].triangleindices_;

    int numBoundary=0;
    for (int i=0;i<3;i++){
        int next = (i+1) % 3;
        INT2 edge = MeshTools::makeEdge(vInd[i],vInd[next]);
        std::vector<int> faces;
        bool found = Algorithm::FindInMap(mapEdgeToFace, edge, faces);
        ASSERT((found), "The edge " << edge << " is not found.");

        if (faces.size() == 1){
            numBoundary += 1;
        }
    }

    if (numBoundary == 2){
        return true;
    }

    return false;
}

void MeshTools::ReconstructMeshAfterFaceCut(Mesh& mesh)
{
    std::vector<std::vector<int>> neighborIndex;
    MeshTools::CalculateVertexNeighbors(mesh, neighborIndex);

    const auto& v = mesh.getvertices();
    const auto& f = mesh.gettriangles();

    int index=0;
    std::vector<int> MapOldToNewIndex(v.size(),-1);
    for (int i=0;i<neighborIndex.size();i++){
        if (neighborIndex[i].size() != 0){
            MapOldToNewIndex[i] = index;            
            index++;
        }
    }

    std::vector<triangle> newFace;
    std::vector<vertex> newVerts;
    for (int i=0;i<f.size();i++){
        triangle t;
        for (int j=0;j<3;j++){
            int oldIndex = f[i].triangleindices_[j];
            int newIndex = MapOldToNewIndex[oldIndex];
            t.triangleindices_[j] = newIndex;        
        }
        newFace.push_back(t);
    }

    for (int i=0;i<v.size();i++){
        if (neighborIndex[i].size() != 0){
            newVerts.push_back(v[i]);
        }
    }

    auto& oldT = mesh.accesstriangles();
    auto& oldV = mesh.accessvertices();
    oldT.clear(); oldT.insert(oldT.end(), newFace.begin(), newFace.end());
    oldV.clear(); oldV.insert(oldV.end(), newVerts.begin(), newVerts.end());
}


void MeshTools::FindTriangleAngles(const Mesh& mesh, std::vector<Real>& angles){
    angles.clear();

    const auto& tri = mesh.gettriangles();
    const auto& vert = mesh.getvertices();

    #pragma omp parallel
    {
        std::vector<Real> localAngles;
        #pragma omp for
        for (int i=0;i<tri.size();i++){
            triangle ti = tri[i];
            // for each triangle we have 3 angles 
            for (int j=0;j<3;j++){
                int next = (j+1) % 3;
                int before  = j-1;

                if (before < 0){
                    before  += 3;
                }

                // find the 2 vectors of the triangle
                Real3 vec1, vec2;
                Real dist1, dist2;
                mesh.getVertexDistance(vert[ti[next]], vert[ti[j]], vec1, dist1);
                mesh.getVertexDistance(vert[ti[before]], vert[ti[j]], vec2, dist2);
                
                Real angle = LinAlg3x3::findAngle(vec1,vec2) * 180 / Constants::PI;
                localAngles.push_back(angle);
            }
        }

        #pragma omp critical
        {
            angles.insert(angles.end(), localAngles.begin(), localAngles.end());
        }
    }
}

void MeshTools::FindSideLengths(const Mesh& mesh, std::vector<Real>& SideLengths){
    const auto& tri = mesh.gettriangles();
    const auto& v = mesh.getvertices();

    SideLengths.clear();

    #pragma omp parallel
    {
        std::vector<Real> localSideLength;

        #pragma omp for       
        for (int i =0;i<tri.size();i++){
            auto t = tri[i];

            // first side length 0 - 1
            Real side1, side2, side3;
            Real3 vec1, vec2, vec3;

            mesh.getVertexDistance(v[0], v[1], vec1, side1);
            mesh.getVertexDistance(v[1], v[2], vec2, side2);
            mesh.getVertexDistance(v[2], v[0], vec3, side3);

            localSideLength.push_back(side1);
            localSideLength.push_back(side2);
            localSideLength.push_back(side3);
        }

        #pragma omp critical
        {
            SideLengths.insert(SideLengths.end(), localSideLength.begin(), localSideLength.end());
        }
    }
}

std::vector<int> MeshTools::RemoveIsolatedVertices(Mesh& mesh){
    std::vector<std::vector<int>> MapVerticesToFaces, neighborIndices;
    MeshTools::MapVerticesToFaces(mesh, MapVerticesToFaces);
    MeshTools::CalculateVertexNeighbors(mesh, neighborIndices);

    std::vector<vertex> newV;
    std::vector<triangle> newF;
    const auto& cv = mesh.getvertices();
    const auto& cf = mesh.gettriangles();
    std::vector<int> MapOldIndexToNew(cv.size(), -100);
    int count = 0;
    for (int i=0;i<MapVerticesToFaces.size();i++){
        // check if this is an isolated triangle --> for all of its neighbors do they map to 1 triangle only
        if (MapVerticesToFaces[i].size() != 0){
            newV.push_back(cv[i]);
            MapOldIndexToNew[i] = count;
            count ++;
        }
    }

    for (int i=0;i<cf.size();i++){
        triangle t;
        for (int j=0;j<3;j++){
            t[j] = MapOldIndexToNew[cf[i][j]];
            ASSERT((t[j] != -100), "Something went wrong.");
        }
        if (t[0] == t[1] || t[1] == t[2] || t[2] == t[0]){
            std::cout << "Something went wrong in Removing isolated vertice " << std::endl;
        }
        newF.push_back(t);
    }

    auto& v = mesh.accessvertices();
    auto& f = mesh.accesstriangles();
    v.clear();
    v.insert(v.end(), newV.begin(), newV.end());
    f.clear();
    f.insert(f.end(), newF.begin(), newF.end());

    return MapOldIndexToNew;
}

void MeshTools::RemoveIsolatedFaces(Mesh& mesh){
    std::map<INT2, std::vector<int>> EtoF;
    std::vector<std::vector<INT2>> VtoE;
    MeshTools::MapEdgeToFace(mesh, EtoF, VtoE, false);

    std::vector<triangle> newF;
    const auto& cf = mesh.gettriangles();
    for (int i=0;i<cf.size();i++){
        if (! IsIsolatedFace(mesh, i, EtoF)){
            newF.push_back(cf[i]);
        }
    }

    auto& f = mesh.accesstriangles();
    f.clear();
    f.insert(f.end(), newF.begin(), newF.end());
}

void MeshTools::RemoveDuplicatedFaces(Mesh& mesh){
    std::map<INT3, std::vector<int>> MapVertexIndexToFace;

    const auto& f = mesh.gettriangles();
    const auto& v = mesh.getvertices(); 

    for (int i=0;i<f.size();i++){
        INT3 indices = f[i].triangleindices_;
        // sort the indices
        Algorithm::sort(indices);
        Algorithm::InsertInVectorMap(MapVertexIndexToFace, indices,i);
    }

    std::vector<triangle> newf;
    for (auto it = MapVertexIndexToFace.begin(); it != MapVertexIndexToFace.end(); it++){
        // unique triangle
        if (it -> second.size() == 1){
            int t_index = it -> second[0];
            newf.push_back(f[t_index]);
        }
        else{
            std::vector<int> tempVec = it -> second;
            Algorithm::sort(tempVec);
            int t_index = tempVec[0];
            newf.push_back(f[t_index]);
        }
    }

    auto& nf = mesh.accesstriangles();
    nf.clear();
    nf.insert(nf.end(), newf.begin(), newf.end());
}

void MeshTools::RemoveMinimumNeighbors(Mesh& m, int num_search, int min_num_neighbors)
{
    std::vector<std::vector<int>> vertex_neighbors;
    MeshTools::CalculateVertexNeighbors(m, vertex_neighbors);

    std::vector<std::vector<int>> NearbyNeighbors;
    Graph::BFS_kring_neighbor(vertex_neighbors, num_search, NearbyNeighbors);

    // check whether or not the nearby neighbors exceed a certain threshold
    const auto& Oldv = m.getvertices();
    const auto& Oldt = m.gettriangles();
    std::vector<int> MapOldToNew(Oldv.size(),-1);
    int index = 0;

    // construct a mapping from old to new
    for (int i=0;i<NearbyNeighbors.size();i++){
        if (NearbyNeighbors[i].size() >= min_num_neighbors){
            MapOldToNew[i] = index;
            index++;
        }
    }

    std::vector<vertex> Newv;
    std::vector<triangle> Newt;
    // then let construct the new vertex and triangles
    for (int i=0;i<Oldv.size();i++){
        if (MapOldToNew[i] != -1){
            Newv.push_back(Oldv[i]);
        }
    }

    for (int i=0;i<Oldt.size();i++){
        bool isValidTriangle=true;
        triangle newtriangle;
        for (int j=0;j<3;j++){
            if (MapOldToNew[Oldt[i][j]] == -1){
                isValidTriangle = false;
            }
            else{
                newtriangle[j] = MapOldToNew[Oldt[i][j]];
            }
        }

        if (isValidTriangle){
            Newt.push_back(newtriangle);
        }
    }

    m.SetVerticesAndTriangles(Newv, Newt);
}

void MeshTools::ChangeWindingOrderPosZ(Mesh& m){
    std::vector<Real> Area;
    std::vector<Real3> FN;
    MeshTools::CalculateTriangleAreasAndFaceNormals(m, Area, FN);

    for (int i=0;i<FN.size();i++){
        if (FN[i][2] < 0){
            ChangeWindingOrder(m, i);
        }
    }
}

void MeshTools::ChangeWindingOrder(Mesh& m){
    for (int i=0;i<m.gettriangles().size();i++){
        ChangeWindingOrder(m, i);
    }
}


void MeshTools::ChangeWindingOrder(Mesh& m, int num){
    auto& tri = m.accesstriangles();

    triangle nT;
    nT[0] = tri[num][1];
    nT[1] = tri[num][0];
    nT[2] = tri[num][2];

    tri[num] = nT;
}

void MeshTools::MeshPlaneClipping(Mesh& m, Real3& point, Real3& plane){
    // first we need to find the signed distance of all the vertex points to the plane 
    const auto& ori_v = m.getvertices();
    const auto& ori_f = m.gettriangles();
    int num_v = ori_v.size();
    int num_f = ori_f.size();
    std::vector<Real> signed_distance(num_v, 0.0);
    std::vector<int> vertex_lower;
    for (int i=0;i<num_v;i++){
        Real3 diff = ori_v[i].position_ - point; 
        Real dist  = LinAlg3x3::DotProduct(diff, plane);
        signed_distance[i] = dist;
        if (dist < 0){
            vertex_lower.push_back(i);
        }
    }

    // then figure out which triangles are intersected by the plane 
    std::vector<int> intersect_triangle;
    std::vector<int> kept_triangle;
    std::vector<std::vector<int>> edges_intersected;
    std::vector<uint8_t> venums(num_f);

    for (auto&& [index, f] : Algorithm::enumerate(ori_f)){
        Real dist1 = signed_distance[f[0]]; Real dist2 = signed_distance[f[1]]; Real dist3 = signed_distance[f[2]];
        uint8_t mask;

        if ((dist1 < 0) && (dist2 < 0) && (dist3<0)){
            kept_triangle.push_back(index);
        }
        else if ((dist1 > 0) && (dist2 > 0) && (dist3 >0)){
            kept_triangle.push_back(index);
        }
        else{
            intersect_triangle.push_back(index);
        }

        // write to bits
        if (dist1 < 0) mask |= 1 << 0;
        if (dist2 < 0) mask |= 1 << 1;
        if (dist3 < 0) mask |= 1 << 2;

        venums[index] = mask;
    }

    // figure out how many total triangles there are
    std::array<int,8> s_clipVertCountTable{0,1,1,2,1,2,2,1};
    int X = std::numeric_limits<int>::max();
    std::array<std::array<int,6>,8> s_clipTriTable{{
        {X,X,X,X,X,X}, 
        {0,3,5,X,X,X}, 
        {3,1,4,X,X,X}, 
        {0,1,5,1,4,5}, 
        {4,2,5,X,X,X}, 
        {0,3,4,0,4,2}, 
        {1,5,3,1,2,5},
        {0,1,2,X,X,X}
    }};

    int new_num_f= 0;
    for (int i=0;i<num_f;i++){
        new_num_f += s_clipVertCountTable[venums[i]];
    }

}

void MeshTools::TriangleCases(std::vector<INT3>& signs, std::vector<bool>& basic, std::vector<bool>& one_vertex, std::vector<bool>& one_edge){
    /*
    code : signs      : intersects
    0    : [-1 -1 -1] : No
    2    : [-1 -1  0] : No
    4    : [-1 -1  1] : Yes; 2 on one side, 1 on the other
    6    : [-1  0  0] : Yes; one edge fully on plane
    8    : [-1  0  1] : Yes; one vertex on plane 2 on different sides
    12   : [-1  1  1] : Yes; 2 on one side, 1 on the other
    14   : [0 0 0]    : No (on plane fully)
    16   : [0 0 1]    : Yes; one edge fully on plane
    20   : [0 1 1]    : No
    28   : [1 1 1]    : No
    */
    std::vector<int> coded(signs.size(),14);
    std::vector<bool> key(29, false);
    one_edge.clear(); one_edge.resize(signs.size(), false);
    one_vertex.clear(); one_vertex.resize(signs.size(), false);
    basic.clear(); basic.resize(signs.size(), false);

    #pragma omp parallel for
    for (int i=0;i<signs.size();i++){
        Algorithm::sort(signs[i]);
    }

    #pragma omp parallel for 
    for (int i=0;i<signs.size();i++){
        for (int j=0;j<3;j++){
            coded[i] += (signs[i][j] << (3-j));
        }
    }

    #pragma omp parallel for
    for (int i=0;i<signs.size();i++){
        if (coded[i] == 16){
            one_edge[i] = true;
        }
        else if (coded[i] == 8){
            one_vertex[i] = true;
        }
        else if ((coded[i] == 4) || (coded[i]=12)){
            basic[i] = true;
        }
    }
}

void MeshTools::MeshPlaneIntersection(Mesh& m, Real3& points, Real3& normal){
    const auto& v = m.getvertices();
    const auto& f = m.gettriangles();
    std::vector<Real> dots(v.size(), 0.0);
    std::vector<int> signs(v.size(), 0);
    std::vector<INT3> face_signs(f.size());
    Real tolerance = 1e-5;

    #pragma omp parallel for
    for (int i=0;i<v.size();i++){
        dots[i] = LinAlg3x3::DotProduct(v[i].position_ - points, normal);
        if (dots[i] > tolerance){
            signs[i] = 1;
        }
        else if (dots[i] < - tolerance){
            signs[i] = -1;
        }
        else{
            signs[i] = 0;
        }
    }

    // get the face signs ready 
    #pragma omp parallel for
    for (int i=0;i<f.size();i++){
        for (int j=0;j<3;j++){
            int index = f[i][j];
            face_signs[i][j] = signs[index];
        }
    }

    // get basic, one_vertex and one_edge
    std::vector<bool> basic, one_vertex, one_edge;
    TriangleCases(face_signs, basic, one_vertex, one_edge);
}

void MeshTools::CalculateCotangentWeights(Mesh& m, const std::vector<std::vector<int>>& neighborIndices, const std::map<INT2, std::vector<int>>& MapEdgeToFace, const std::map<INT2, std::vector<int>>& MapEdgeToOpposingVerts, std::vector<Real3>& dAdpi)
{
    // obtain triangles and vertices from mesh
    const auto& triangles = m.gettriangles();
    const auto& vertices  = m.getvertices();
    dAdpi.clear();
    dAdpi.resize(vertices.size(), {});

    // get the neighbor indices 
    #pragma omp parallel for
    for (int i=0;i<neighborIndices.size();i++){
        // find the neighbor of vertex i
        std::vector<int> neighbors = neighborIndices[i];

        // get neighbor size 
        int neighborSize = neighbors.size();

        Real3 dAdpi_j ={0,0,0};

        // iterate over the neighbor indices
        for (int j=0;j<neighborSize;j++){
            // fi - fj
            Real3 vec_fj_fi;
            Real vec_fj_fi_sq;

            // fi - fj
            m.getVertexDistance(vertices[i], vertices[neighbors[j]],  vec_fj_fi, vec_fj_fi_sq);

            // make an edge between this vertex and its neighbor
            INT2 edge = MeshTools::makeEdge(i,neighbors[j]);

            // map this edge to the faces that it corresponds to 
            std::vector<int> faces;
            bool FoundEdge = Algorithm::FindInMap(MapEdgeToFace, edge, faces);
            ASSERT(FoundEdge, "The edge " << edge <<  " is not found.");

            // map this edge to the opposing vertices indices 
            std::vector<int> OpposingVerts;
            bool FoundVerts = Algorithm::FindInMap(MapEdgeToOpposingVerts, edge, OpposingVerts);
            ASSERT(FoundVerts, "The edge " << edge << " is not found.");

            // factor here is cot(\alpha) + cot(\beta) --> where \alpha and \beta are the opposing angles shared by an edge 
            Real cotOpposingAnglesSum = 0.0;
            int index1 = edge[0];
            int index2 = edge[1];

            // calculate the cosine theta with respect to opposing points 
            for (int OpposingIdx : OpposingVerts){
                Real3 vec1, vec2;
                Real vec1sq, vec2sq;

                // calculate the vertex distance and vector 
                m.getVertexDistance(vertices[OpposingIdx], vertices[index1], vec1, vec1sq);
                m.getVertexDistance(vertices[OpposingIdx], vertices[index2], vec2, vec2sq);

                // normalize the vector
                LinAlg3x3::normalize(vec1);
                LinAlg3x3::normalize(vec2);

                // find the costine angle between
                Real costheta  = LinAlg3x3::findCosangle(vec1, vec2);
                Real sin2theta = 1 - costheta*costheta;

                if (sin2theta < 1e-8){
                    std::cout << "sin2theta = " << sin2theta << std::endl;
                    std::cout << "triangle " << index1 << " " << index2 << " " << OpposingIdx << " is wrong" << std::endl;
                }

                Real sintheta = std::sqrt(sin2theta);
                cotOpposingAnglesSum += costheta/sintheta;
            }

            // add to 
            dAdpi_j = dAdpi_j + cotOpposingAnglesSum * vec_fj_fi; 
        }

        // add to dAdpi
        dAdpi[i] = 0.5 * dAdpi_j;
    }
}

void MeshTools::CalculateAreaDerivatives(Mesh& m, std::vector<Real3>& dAdpi){
    std::vector<std::vector<int>> neighborIndices;
    std::map<INT2, std::vector<int>> etof;
    std::map<INT2, std::vector<int>> etoo;

    CalculateVertexNeighbors(m, neighborIndices);
    MapEdgeToFace(m, etof);
    MapEdgeToOpposingVertices(m, etof, etoo);
    CalculateCotangentWeights(m, neighborIndices, etof, etoo, dAdpi);
}


void MeshTools::CalculateVolumeDerivatives(Mesh& m, const std::vector<std::vector<int>>& vtof, std::vector<Real3>& VolumeDerivatives, Real3 shift){
    const auto& verts = m.getvertices();
    const auto& faces = m.gettriangles();

    VolumeDerivatives.clear();
    VolumeDerivatives.resize(verts.size());

    #pragma omp parallel for
    for (int i=0;i<vtof.size();i++){
        // initialize the gradient 
        Real3 gradient = {0,0,0};

        // iterate over all face neighbors
        for (int face_idx : vtof[i]){
            // find the vertices that is i
            int triangle_indices;

            // iterate over the triangle indices to find it
            for (int j=0;j<3;j++){
                if (faces[face_idx].triangleindices_[j] == i){
                    triangle_indices = j;
                    break;
                }
            }

            // find the other 2 vertex indices 
            int idx1,idx2;
            Real3 shifted_p1, shifted_p2;
            idx1 = faces[face_idx][(triangle_indices + 1)%3];
            idx2 = faces[face_idx][(triangle_indices + 2)%3];

            // find the shifted position of the other 2 vertices with respect to the first vertex 
            shifted_p1 = m.getShiftedVertexPosition(verts[idx1], verts[i]) + shift;
            shifted_p2 = m.getShiftedVertexPosition(verts[idx2], verts[i]) + shift;

            // find the gradient
            gradient = gradient + 1.0 / 6.0 * LinAlg3x3::CrossProduct(shifted_p1, shifted_p2);
        }

        VolumeDerivatives[i] = gradient;
    }
}

void MeshTools::CalculateVolumeDerivatives(Mesh& m, std::vector<Real3>& VolumeDerivatives, Real3 shift){
    std::vector<std::vector<int>> vtof;
    MapVerticesToFaces(m, vtof);

    CalculateVolumeDerivatives(m, vtof, VolumeDerivatives, shift);
}

MeshTools::Real MeshTools::CalculateVolumeDivergenceTheorem(Mesh& m, const std::vector<Real>& vecArea, const std::vector<Real3>& Normal){
    // use divergence theorem to calculate the volume of a mesh
    const auto& verts = m.getvertices();
    const auto& faces = m.gettriangles();

    Real volume = 0.0f;

    // iterate over faces
    for (int i=0;i<faces.size();i++){
        // first calculate the center
        Real3 diff1, diff2, centroid;
        Real diffsq1, diffsq2;
        Real3 v0,v1,v2;

        // initialize the positions
        v0 = verts[faces[i][0]].position_; v1 = verts[faces[i][1]].position_; v2 = verts[faces[i][2]].position_;
        

        // shift vertices 2 and 3 wrt to the first 
        if (m.isPeriodic()){
            m.getVertexDistance(verts[faces[i][1]], verts[faces[i][0]], diff1, diffsq1);
            m.getVertexDistance(verts[faces[i][2]], verts[faces[i][0]], diff2, diffsq2);

            // calculate the centroid
            centroid = 1.0 / 3.0 * (v0 + v0 + diff1 + v0 + diff2);

            // shift centroid to center of box
            Real3 shift = m.getShiftIntoBox(centroid);
            centroid = centroid + shift;
        }
        else{
            centroid = 1.0 / 3.0 * (v0 + v1 + v2);
        }

        // get the volume
        volume += 1.0 / 3.0 * LinAlg3x3::DotProduct(centroid, Normal[i]) * vecArea[i];
    }

    return volume;
}


MeshTools::Real MeshTools::CalculateArea(Mesh& m, std::vector<Real>& vecArea, std::vector<Real3>& normals){
    // normals and ares
    vecArea.clear();
    normals.clear();

    // calculate the triangle area and face normals 
    MeshTools::CalculateTriangleAreasAndFaceNormals(m, vecArea, normals);

    Real sum_area=0.0;

    for (int i=0;i<vecArea.size();i++){
        sum_area += vecArea[i];
    }

    return sum_area;
}

void MeshTools::CVT_optimize_Mesh(Mesh& m, Real volume_weight, int nb_iterations){
    std::vector<Real3> v = m.getVertexPositions();
    std::vector<INT3>  f = m.getFaces();

    // start a edge mesh
    CMC::EdgeMesh emesh(v,f);
    emesh.mark_boundary();
    CMC::CMC_Evolver evolver;
    evolver.set_mesh(&emesh);
    evolver.set_volume_weight(volume_weight);
    evolver.cmc_qnewton(nb_iterations,1);

    std::vector<Real3> new_v;
    std::vector<INT3> new_f;
    emesh.get_verts_faces(new_v, new_f);
    Mesh new_m(new_v, new_f);

    m = new_m;
}

void MeshTools::CGAL_optimize_Mesh(Mesh& m, int nb_iterations, Real degree, bool use_restriction, bool AddNewTriangles){
    M cgal_m;
    Real3 box;
    bool isPeriodic = m.isPeriodic();

    if (isPeriodic){
        MeshTools::ConvertToNonPBCMesh(m, AddNewTriangles);
        box = m.getBoxLength();
    }

    std::vector<std::array<FT,3>> points;
    std::vector<CGAL_Polygon> polygons;

    const auto& verts = m.getvertices();
    const auto& face = m.gettriangles();

    for (auto v : verts){
        points.push_back(CGAL::make_array<FT>(v.position_[0], v.position_[1], v.position_[2]));
    }

    for (auto f : face){
        polygons.push_back({(std::size_t)f.triangleindices_[0], (std::size_t)f.triangleindices_[1], (std::size_t)f.triangleindices_[2]});
    }

    // repair polygon 
    PMP::repair_polygon_soup(points, polygons, CGAL::parameters::geom_traits(Array_traits()));
    PMP::orient_polygon_soup(points, polygons);
    PMP::polygon_soup_to_polygon_mesh(points, polygons, cgal_m);

    typedef boost::property_map<M, CGAL::edge_is_feature_t>::type EIFMap;
    EIFMap eif = get(CGAL::edge_is_feature, cgal_m);
    PMP::detect_sharp_edges(cgal_m, degree, eif);

    PMP::angle_and_area_smoothing(cgal_m, CGAL::parameters::number_of_iterations(nb_iterations)
                                                       .use_safety_constraints(use_restriction) // authorize all moves
                                                       .edge_is_constrained_map(eif));

    std::vector<Real3> new_p_list; 
    std::vector<INT3>  new_f_list;
    for(auto v : cgal_m.vertices()){
        auto p = cgal_m.point(v);
        Real3 new_p = {p[0], p[1], p[2]};
        new_p_list.push_back(new_p);
    }

    for (auto f : cgal_m.faces()){
        std::vector<int> new_f_vec;
        for (auto v : CGAL::vertices_around_face(cgal_m.halfedge(f), cgal_m)){
            new_f_vec.push_back(v.id());
        }
        INT3 new_f;
        new_f[0] = new_f_vec[0]; new_f[1] = new_f_vec[1]; new_f[2] = new_f_vec[2];
        new_f_list.push_back(new_f);
    }

    Mesh new_m(new_p_list, new_f_list);

    if (isPeriodic){
        new_m.setBoxLength(box);
        MeshTools::MakePBCMesh(new_m);
    }

    m = new_m;
    m.CalcVertexNormals();
}

void MeshTools::CalculateAVnbs(Mesh& m, AFP_shape* s, Real& A, Real& V, int v_num, bool useNumerical, Real3 Vshift){
    std::vector<int> BoundaryIndices;
    std::vector<Real> ulist, vlist;

    // obtain boundary indices 
    MeshTools::CalculateBoundaryVerticesIndex(m, BoundaryIndices);

    // find the ulist and vlist --> this also rearranges boundary indices 
    MeshTools::FindBoundaryUV(m, ulist, vlist, BoundaryIndices, s,true);
    MeshTools::CalculateAVnbs(m, s, BoundaryIndices, ulist, vlist, A, V, v_num, useNumerical, Vshift);
}


void MeshTools::CalculateAVnbs(Mesh& m, AFP_shape* s,std::vector<int>& BoundaryIndices,\
                               std::vector<Real>& ulist, std::vector<Real>& vlist, Real& A, Real& V,\
                               int num_v, bool useNumerical, Real3 Vshift){
    // We assume that the ulist is ordered
    A = 0.0;
    V = 0.0;
    int numBoundary = ulist.size();
    const auto& verts = m.getvertices();
    for (int i=0;i<vlist.size();i++){
        Real vstep = (Constants::PI/2 - vlist[i]) / num_v;
        Real u = ulist[i];
        Real du = ulist[(i+1) % numBoundary] - ulist[i];
        if (du < 0){
            du += 2 * Constants::PI;
        }

        int ind = BoundaryIndices[i];

        #pragma omp parallel 
        {
            Real A_local= 0.0;
            Real V_local= 0.0;
            #pragma omp for
            for (int j=0;j<num_v;j++){
                Real3 drdu, drdv;
                Real v = vlist[i] + j * vstep;

                drdu = s->drdu(u,v, useNumerical);
                drdv = s->drdv(u,v, useNumerical);

                Real3 cp = LinAlg3x3::CrossProduct(drdu, drdv);
                
                A_local += std::sqrt(LinAlg3x3::DotProduct(cp,cp)) * du * vstep;

                Real3 pos = verts[ind].position_ + Vshift;

                V_local += 1.0/3.0 * LinAlg3x3::DotProduct(pos, cp) * du * vstep;
            }

            #pragma omp critical
            {
                A += A_local;
                V += V_local;
            }
        }
    }
}

std::unique_ptr<AFP_shape> MeshTools::ReadAFPShape(CommandLineArguments& cmd){
    std::string shape_name;
    cmd.readString("shape", CommandLineArguments::Keys::Required, shape_name);
    ParameterPack shapePack;

    // initialize the shape
    std::unique_ptr<AFP_shape> shape;

    if (shape_name == "SuperEgg"){
        std::string a,b,n,zmax,a_taper,b_taper,a_alpha,b_alpha;
        std::vector<std::string> center;
        bool read;
        cmd.readValue("a", CommandLineArguments::Keys::Required, a);
        shapePack.insert("a", a);
        cmd.readValue("b", CommandLineArguments::Keys::Required, b);
        shapePack.insert("b", b);
        cmd.readValue("zmax", CommandLineArguments::Keys::Required, zmax);
        shapePack.insert("zmax", zmax);
        cmd.readVector("center", CommandLineArguments::Keys::Required, center);
        shapePack.insert("center", center);
        read = cmd.readValue("n", CommandLineArguments::Keys::Optional, n);
        if (read){
            shapePack.insert("n", n);
        }
        read = cmd.readValue("a_taper", CommandLineArguments::Keys::Optional, a_taper);
        if (read){
            shapePack.insert("a_taper", a_taper);
        }
        read = cmd.readValue("b_taper", CommandLineArguments::Keys::Optional, b_taper);
        if (read){
            shapePack.insert("b_taper", b_taper);
        }
        read = cmd.readValue("a_alpha", CommandLineArguments::Keys::Optional, a_alpha);
        if (read){
            shapePack.insert("a_alpha", a_alpha);
        }
        read = cmd.readValue("b_alpha", CommandLineArguments::Keys::Optional, b_alpha);
        if (read){
            shapePack.insert("b_alpha", b_alpha);
        }

        // initialize the shape 
        shape = std::make_unique<SuperEgg>(shapePack);
    }
    else if (shape_name == "Sphere"){
        std::string radius;
        std::vector<std::string> center;

        cmd.readVector("center", CommandLineArguments::Keys::Required, center);
        cmd.readValue("radius", CommandLineArguments::Keys::Required, radius);

        shapePack.insert("radius", radius);
        shapePack.insert("center", center);

        shape = std::make_unique<Sphere>(shapePack);
    }
    else if (shape_name == "BulgingSphere"){
        std::string radius;
        std::string bulge_factor;
        std::string theta_phi_factor;
        std::vector<std::string> center;

        cmd.readVector("center", CommandLineArguments::Keys::Required, center);
        cmd.readValue("radius", CommandLineArguments::Keys::Required, radius);
        cmd.readValue("bulge_factor", CommandLineArguments::Keys::Optional, bulge_factor);
        cmd.readValue("theta_phi_factor", CommandLineArguments::Keys::Optional, theta_phi_factor);

        shapePack.insert("radius", radius);
        shapePack.insert("center", center);
        shapePack.insert("bulge_factor", bulge_factor);
        shapePack.insert("theta_phi_factor", theta_phi_factor);

        shape = std::make_unique<BulgingSphere>(shapePack);
    }

    return std::move(shape);
}

void MeshTools::CalculateBoundaryVerticesIndex(const Mesh& mesh, std::vector<int>& boundaryIndices, bool order, Real3 center){
    std::vector<bool> boundaryIndicator;
    boundaryIndices.clear();
    CalculateBoundaryVertices(mesh, boundaryIndicator);

    for (int i=0;i<mesh.getvertices().size();i++){
        if (IsBoundary(i, boundaryIndicator)){
            boundaryIndices.push_back(i);
        }
    }

    if (order){
        const auto& vertices = mesh.getvertices();
        std::vector<Real> theta;
        // shift all the vertices wrt to the first vertex 
        int first_ind = boundaryIndices[0];
        std::vector<Real3> shifted_pos;
        Real3 shifted_pos_sum = {0,0,0};

        for (int ind : boundaryIndices){
            // find the shifted vertex position
            Real3 s_p = mesh.getShiftedVertexPosition(vertices[ind].position_, center);
            shifted_pos.push_back(s_p);
        }

        // find the center
        Real3 center = shifted_pos_sum / (Real)boundaryIndices.size();

        // shift all the position by the center and calculate the atan2
        for (int i=0;i<shifted_pos.size();i++){
            Real t = std::atan2(shifted_pos[i][1], shifted_pos[i][0]);
            if (t < 0.0f){
                t += Constants::PI * 2.0f;
            }
            theta.push_back(t);
        }

        std::vector<int> index = Algorithm::argsort(theta);

        // rewrite the boundary indices 
        std::vector<int> newBoundaryIndex;
        for (int i=0;i<index.size();i++){
            newBoundaryIndex.push_back(boundaryIndices[index[i]]);
        }

        boundaryIndices = newBoundaryIndex;
    }
}


MeshTools::refineptr MeshTools::ReadInterfacialMin(CommandLineArguments& cmd){
    std::string stepsize, ca_file_output="ca.out", T="298", L="0", maxstep="1e5", tolerance="0.00001", printevery="1000", optimizeevery="1e10";
    std::string MaxStepCriteria="true";

    cmd.readString("maxstep", CommandLineArguments::Keys::Optional, maxstep);
    cmd.readString("temperature", CommandLineArguments::Keys::Required, T);
    cmd.readString("stepsize", CommandLineArguments::Keys::Optional, stepsize);
    cmd.readString("tolerance", CommandLineArguments::Keys::Optional, tolerance);
    cmd.readValue("printevery", CommandLineArguments::Keys::Optional, printevery);
    cmd.readString("Lagrange", CommandLineArguments::Keys::Optional, L);
    cmd.readString("optimize_every", CommandLineArguments::Keys::Optional, optimizeevery);

    // define parameter packs
    ParameterPack refinePack;
    refinePack.insert("name", "refine");
    refinePack.insert("optimize_every", optimizeevery);
    refinePack.insert("maxstep", maxstep);
    refinePack.insert("temperature", T);
    refinePack.insert("L", L);
    refinePack.insert("stepsize", stepsize);
    refinePack.insert("tolerance", tolerance);
    refinePack.insert("print_every", printevery);
    refinePack.insert("MaxStepCriteria", MaxStepCriteria);
    
    // initialize refine ptr
    refineptr r;
    MeshRefineStrategyInput input = {{refinePack}};
    r = refineptr(MeshRefineStrategyFactory::Factory::instance().create("InterfacialFE_minimization", input));

    return std::move(r);
}

MeshTools::refineptr MeshTools::ReadInterfacialMinBoundary(CommandLineArguments& cmd){
    std::string stepsize, ca_file_output="ca.out", T="298", L="0", maxstep="1e5", tolerance="0.00001", printevery="1000", optimize_every="1e10";
    std::string boundarymaxstep="10000", boundarystepsize="0.5", boundarytolerance="5e-6", L2tolerance="5e-5", boundary_optimize_every="300";
    std::string MaxStepCriteria="true", boundaryMaxStepCriteria="true", L2="0", L2_step_size="10.0",L2maxstep="1e5", zstar, zstar_deviation="0.003";
    std::string dgamma_gamma, useNumerical="false";

    cmd.readString("maxstep", CommandLineArguments::Keys::Optional, maxstep);
    cmd.readString("temperature", CommandLineArguments::Keys::Required, T);
    cmd.readString("stepsize", CommandLineArguments::Keys::Optional, stepsize);
    cmd.readString("tolerance", CommandLineArguments::Keys::Optional, tolerance);
    cmd.readString("printevery", CommandLineArguments::Keys::Optional, printevery);
    cmd.readString("Lagrange", CommandLineArguments::Keys::Optional, L);
    cmd.readString("optimize_every", CommandLineArguments::Keys::Optional, optimize_every);

    // boundary terms
    cmd.readString("boundarymaxstep", CommandLineArguments::Keys::Optional, boundarymaxstep);
    cmd.readString("boundarystepsize", CommandLineArguments::Keys::Optional, boundarystepsize);
    cmd.readString("boundarytolerance", CommandLineArguments::Keys::Optional, boundarytolerance);
    cmd.readString("boundary_optimize_every", CommandLineArguments::Keys::Optional, boundary_optimize_every);
    cmd.readString("dgamma_gamma", CommandLineArguments::Keys::Required, dgamma_gamma);
    cmd.readString("useNumerical", CommandLineArguments::Keys::Optional, useNumerical);
    cmd.readString("L2_guess", CommandLineArguments::Keys::Optional, L2);
    cmd.readString("L2_step_size", CommandLineArguments::Keys::Optional, L2_step_size);
    cmd.readString("L2_tolerance", CommandLineArguments::Keys::Optional, L2tolerance);
    cmd.readString("L2maxstep", CommandLineArguments::Keys::Optional, L2maxstep);
    cmd.readString("zstar_deviation", CommandLineArguments::Keys::Optional, zstar_deviation);
    cmd.readString("zstar", CommandLineArguments::Keys::Optional, zstar);
    cmd.readString("MaxStepCriteria", CommandLineArguments::Keys::Optional, MaxStepCriteria);
    cmd.readString("boundaryMaxStepCriteria", CommandLineArguments::Keys::Optional, boundaryMaxStepCriteria);

    // define parameter packs
    ParameterPack refinePack;
    refinePack.insert("name", "refine");
    refinePack.insert("optimize_every", optimize_every);
    refinePack.insert("maxstep", maxstep);
    refinePack.insert("temperature", T);
    refinePack.insert("L", L);
    refinePack.insert("stepsize", stepsize);
    refinePack.insert("tolerance", tolerance);
    refinePack.insert("print_every", printevery);
    refinePack.insert("MaxStepCriteria", MaxStepCriteria);
    refinePack.insert("boundaryMaxStepCriteria", boundaryMaxStepCriteria);

    // boundary terms
    refinePack.insert("boundarymaxstep", boundarymaxstep);
    refinePack.insert("boundarystepsize", boundarystepsize);
    refinePack.insert("boundarytolerance", boundarytolerance);
    refinePack.insert("boundary_optimize_every", boundary_optimize_every);
    refinePack.insert("dgamma_gamma", dgamma_gamma);
    refinePack.insert("L2_guess", L2);
    refinePack.insert("zstar", zstar);
    refinePack.insert("zstar_deviation", zstar_deviation);
    refinePack.insert("L2_step_size", L2_step_size);
    refinePack.insert("L2tolerance", L2tolerance);
    
    // initialize refine ptr
    refineptr r;
    MeshRefineStrategyInput input = {{refinePack}};
    r = refineptr(MeshRefineStrategyFactory::Factory::instance().create("InterfacialFE_minimization", input));

    return std::move(r);
}

void MeshTools::MakePBCMesh(Mesh& m){
    // calculate vertices distances between i,j
    const auto& v = m.getvertices();
    const auto& face = m.getFaces();

    // calculate vertex distances
    std::vector<std::vector<int>> overlapping_index;
    overlapping_index.resize(v.size());

    Real threshold = 1e-8;

    // find distances between each of the vertices i,j --> make a map of i,j such that i maps to a vector of j which are close enough to it
    #pragma omp parallel for
    for (int i=0;i<v.size();i++){
        for (int j=0;j<v.size();j++){
            if (j != i){
                Real3 dist_vec;
                Real dist;
                m.getVertexDistance(v[i], v[j], dist_vec, dist);

                if (dist < threshold){
                    overlapping_index[i].push_back(j);
                }
            }
        }
    }

    // convert the adajacency matrix to a map
    // we want all the larger indices to map to smaller indices , e.g. 4->1
    std::map<int, int> map;
    for (int i=0;i<overlapping_index.size();i++){
        auto& indices = overlapping_index[i];

        // check if all indices are larger
        bool larger_than=true;
        for (int j=0;j<indices.size();j++){
            if (indices[j] < i){
                larger_than = false;
                break;
            }
        }

        // if it is then we start mapping
        if (larger_than){
            for (int j=0;j<indices.size();j++){
                Algorithm::InsertInMap(map,indices[j],i);
            }
        }
    }


    // create new faces
    std::vector<INT3> newf;
    // use map to start 
    for (auto f : face){
        int ind1,ind2,ind3;

        if(!Algorithm::FindInMap(map, f[0], ind1)){ind1 = f[0];}
        if(!Algorithm::FindInMap(map, f[1], ind2)){ind2 = f[1];}
        if (!Algorithm::FindInMap(map, f[2], ind3)){ind3=f[2];}

        newf.push_back({ind1,ind2,ind3});
        if (ind1 == ind2 || ind2 == ind3 || ind1 == ind3){
            std::cout << "Something went wrong, the triangle was initially " << f[0] << " " << f[1] << " " << f[2] << " now it is " << \
            ind1 << " " << ind2 << " " << ind3 << std::endl;
        }
    }

    const auto& vert = m.getVertexPositions();

    // create new mesh
    Mesh newm = Mesh(vert, newf);

    if (m.isPeriodic()){
        newm.setBoxLength(m.getBoxLength());
    }

    m = newm;

    // remove isolated vertices 
    MeshTools::RemoveIsolatedVertices(m);
}

void MeshTools::ShiftPBCMesh(Mesh& m, Real3& shift){
    Real3 boxLength = m.getBoxLength();

    auto& vertices = m.accessvertices();

    for (auto& v : vertices){
        v.position_ = v.position_ + shift;

        Real3 s = m.getShiftIntoBox(v.position_);
        v.position_ = v.position_ + s;
    }
}

void MeshTools::ConvertToCGALMesh(Mesh& m, M& cgal_m){
    const auto& verts = m.getvertices();
    const auto& faces = m.gettriangles();
    std::vector<CGAL::SM_Vertex_index> cgal_verts;

    for (auto& v : verts){
        cgal_verts.push_back(cgal_m.add_vertex(CGAL::Epick::Point_3(v.position_[0], v.position_[1], v.position_[2])));
    }

    for (auto& t : faces){
        INT3 ind = t.triangleindices_;
        cgal_m.add_face(cgal_verts[ind[0]], cgal_verts[ind[1]], cgal_verts[ind[2]]);
    }
}

void MeshTools::FindBoundaryUV(Mesh& m, std::vector<Real>& ulist, std::vector<Real>& vlist, const std::vector<bool>& BoundaryIndicator, AFP_shape* shape){
    ulist.clear();
    vlist.clear();

    // find the ulist and vlist
    const auto& verts = m.getvertices();
    for (int i=0;i<m.getNumVertices();i++){
        if (MeshTools::IsBoundary(i, BoundaryIndicator)){
            Real u,v;

            if (m.isPeriodic()){
                u = shape->CalculateU(verts[i].position_, m.getBoxLength());
            }
            else{
                u = shape->CalculateU(verts[i].position_);
            }

            v = shape->CalculateV(verts[i].position_);
            ulist.push_back(u); vlist.push_back(v);
        }
    }
}

void MeshTools::FindBoundaryUV(Mesh& m, std::vector<Real>& ulist, std::vector<Real>& vlist, std::vector<int>& BoundaryIndices, AFP_shape* shape, bool order){
    ulist.clear();
    vlist.clear();

    // find the ulist and vlist
    const auto& verts = m.getvertices();
    for (int i=0;i<BoundaryIndices.size();i++){
        int ind = BoundaryIndices[i];

        Real u,v;
        if (m.isPeriodic()){
            u = shape->CalculateU(verts[ind].position_, m.getBoxLength());
        }
        else{
            u = shape->CalculateU(verts[ind].position_);
        }

        v = shape->CalculateV(verts[ind].position_);
        ulist.push_back(u); vlist.push_back(v);
    }

    if (order){
        // let's argsort the ulist
        std::vector<int> order = Algorithm::argsort(ulist);
        std::vector<Real> new_ulist, new_vlist;
        std::vector<int> new_BoundaryIndices;

        for (int i=0;i<order.size();i++){
            new_BoundaryIndices.push_back(BoundaryIndices[order[i]]);
            new_ulist.push_back(ulist[order[i]]);
            new_vlist.push_back(vlist[order[i]]);
        }

        ulist = new_ulist;
        vlist = new_vlist;
        BoundaryIndices = new_BoundaryIndices;
    }
}

void MeshTools::CalculatedrdUV(Mesh& m, AFP_shape* shape, std::vector<int>& BoundaryIndices, std::vector<Real>& ulist, std::vector<Real>& vlist, std::vector<Real3>& drdu, std::vector<Real3>& drdv, bool useNumerical){
    // obtain boundary indices 
    MeshTools::CalculateBoundaryVerticesIndex(m, BoundaryIndices);

    // find the ulist and vlist --> this also rearranges boundary indices 
    MeshTools::FindBoundaryUV(m, ulist, vlist, BoundaryIndices, shape,true);

    drdu.clear(); drdv.clear();

    for (int i=0;i<ulist.size();i++){
        Real3 drdu_this, drdv_this;
        if(useNumerical){
            drdu_this = shape->Numericaldrdu(ulist[i], vlist[i]);
            drdv_this = shape->Numericaldrdv(ulist[i], vlist[i]);
        }
        else{
            drdu_this = shape->Analyticaldrdu(ulist[i], vlist[i]);
            drdv_this = shape->Analyticaldrdv(ulist[i], vlist[i]);
        }

        drdu.push_back(drdu_this);
        drdv.push_back(drdv_this);
    }
}

void MeshTools::CalculatedrdUV(Mesh& m, AFP_shape* shape, std::vector<Real3>& drdu, std::vector<Real3>& drdv, bool useNumerical){
    // obtain boundary indices 
    std::vector<int> BoundaryIndices;
    std::vector<Real> ulist, vlist;
    MeshTools::CalculatedrdUV(m, shape, BoundaryIndices, ulist, vlist, drdu, drdv, useNumerical);
}

void MeshTools::CalculatedAVnbsdUV(Mesh& m, AFP_shape* shape, std::vector<Real2>& dAnbsduv, std::vector<Real2>& dVnbsduv, bool useNumerical, Real3 shift){
    // find the ulist and vlist --> this also rearranges boundary indices 
    std::vector<int> BoundaryIndices;
    std::vector<Real> ulist, vlist;
    std::vector<Real3> drdu, drdv;
    MeshTools::CalculatedAVnbsdUV(m, shape, BoundaryIndices, ulist, vlist, drdu, drdv, dAnbsduv, dVnbsduv, useNumerical, shift);
}

void MeshTools::CalculatedAVnbsdUV(Mesh& m, AFP_shape* shape, std::vector<int>& BoundaryIndices,\
                                   std::vector<Real2>& dAnbsduv, std::vector<Real2>& dVnbsduv,\
                                   bool useNumerical, Real3 shift){
    // find the ulist and vlist --> this also rearranges boundary indices 
    std::vector<Real> ulist, vlist;
    std::vector<Real3> drdu, drdv;
    MeshTools::CalculatedAVnbsdUV(m, shape, BoundaryIndices, ulist, vlist, drdu, drdv, dAnbsduv, dVnbsduv, useNumerical, shift);
}

void MeshTools::CalculatedAVnbsdUV(Mesh& m,AFP_shape* shape, std::vector<int>& BoundaryIndices, \
                                   std::vector<Real3>& drdu, std::vector<Real3>& drdv, std::vector<Real2>& dAnbsduv,\
                                   std::vector<Real2>& dVnbsduv, bool useNumerical, Real3 shift){
    std::vector<Real> ulist, vlist;
    MeshTools::CalculatedAVnbsdUV(m, shape, BoundaryIndices, ulist, vlist, drdu, drdv, dAnbsduv, dVnbsduv, useNumerical, shift);
}


void MeshTools::CalculatedAVnbsdUV(Mesh& m,AFP_shape* shape, std::vector<int>& BoundaryIndices, std::vector<Real>& ulist, \
                                   std::vector<Real>& vlist, std::vector<Real3>& drdu, std::vector<Real3>& drdv, \
                                   std::vector<Real2>& dAnbsduv, std::vector<Real2>& dVnbsduv, bool useNumerical, Real3 shift){
    // find the ulist and vlist --> this also rearranges boundary indices 
    MeshTools::CalculatedrdUV(m, shape, BoundaryIndices, ulist, vlist, drdu, drdv, useNumerical);

    // clear the output
    dAnbsduv.clear();
    dVnbsduv.clear();

    // get the vertices 
    const auto& verts = m.getvertices(); 

    // start iterating over the boundary indices
    for (int i=0;i<BoundaryIndices.size();i++){
        // define int
        int prev_ind, ind, next_ind;

        if (i == 0){
            prev_ind = BoundaryIndices[BoundaryIndices.size()-1];
        }
        else{prev_ind = BoundaryIndices[i-1];}
        ind = BoundaryIndices[i];
        next_ind = BoundaryIndices[(i+1) % BoundaryIndices.size()];

        // calculate distance 
        Real3 diff1, diff2;
        Real dist1,dist2;
        m.getVertexDistance(verts[next_ind], verts[ind], diff1, dist1);
        m.getVertexDistance(verts[ind], verts[prev_ind], diff2, dist2);

        // calculate dAnbsdu and dAnbsdv --> keep drdv the same 
        Real3 crossP1          = 0.5 * LinAlg3x3::CrossProduct(diff1, drdv[i]);
        Real3 crossP2          = 0.5 * LinAlg3x3::CrossProduct(diff2, drdv[i]);
        Real dAnbsdv_this      = - LinAlg3x3::norm(crossP1) - LinAlg3x3::norm(crossP2);
        Real dVnbsdv_this      = -1.0 / 3.0 * LinAlg3x3::DotProduct(verts[ind].position_ + shift, crossP1) \
                                 -1.0 / 3.0 * LinAlg3x3::DotProduct(verts[ind].position_ + shift, crossP2); 
        dAnbsduv.push_back({0.0, dAnbsdv_this});
        dVnbsduv.push_back({0.0, dVnbsdv_this});
    }
}


void MeshTools::CalculateContactAngle(Mesh& m, AFP_shape* shape, std::vector<Real>& ca){
    ca.clear();

    // for the shape, we calculate s --> obtain u and v 
    std::vector<int> BoundaryIndices;
    std::vector<Real> ulist, vlist;
    MeshTools::CalculateBoundaryVerticesIndex(m,  BoundaryIndices);

    // in this function, we also ordered the boundaryindices
    MeshTools::FindBoundaryUV(m, ulist, vlist, BoundaryIndices, shape, true);

    // now we calculate s,t etc.
    const auto& verts = m.getvertices();
    for (int i=0;i<BoundaryIndices.size();i++){
        int ind = BoundaryIndices[i];

        // get the normals
        auto N  = verts[ind].normals_;

        // get the s 
        Real3 s,t;
        shape->CalculateNumericalNormalAndTangent(ulist[i], vlist[i],t,s);

        // find the dot product between s and N
        ca.push_back(LinAlg3x3::DotProduct(N,s));
    }
}

void MeshTools::CalculateContactAngleDerivative(Mesh& m, AFP_shape* shape, \
                                std::vector<Real>& ca_list, Real k0, Real3 Volume_shift, bool useNumerical){
    // calculate area derivative --> dAdr at the boundary and volume derivative dVdr 
    std::vector<Real3> dAdr ,dVdr;
    MeshTools::CalculateAreaDerivatives(m, dAdr);
    MeshTools::CalculateVolumeDerivatives(m, dVdr, Volume_shift);

    // calculate drdu, drdv, boundaryindices, dAnbsdv, dAnbsdu, dVnbsdv, dVnbsdu
    std::vector<int> BoundaryIndices;
    std::vector<Real2> dAnbsduv, dVnbsduv;
    std::vector<Real3> drdu, drdv;
    std::vector<Real> ulist, vlist;
    MeshTools::CalculatedAVnbsdUV(m, shape, BoundaryIndices, ulist, \
                                vlist, drdu, drdv, dAnbsduv, dVnbsduv, useNumerical, Volume_shift);
    ca_list.clear();

    for (int j=0;j<BoundaryIndices.size();j++){
        // get the actual index of boundary
        int ind = BoundaryIndices[j];

        Real dAdu = LinAlg3x3::DotProduct(drdu[j], dAdr[ind]);
        Real dAdv = LinAlg3x3::DotProduct(drdv[j], dAdr[ind]);

        Real dVdu = LinAlg3x3::DotProduct(drdu[j], dVdr[ind]);
        Real dVdv = LinAlg3x3::DotProduct(drdv[j], dVdr[ind]);

        // calculate dAnbsdu and dAnbsdv --> keep drdv the same 
        Real dAnbsdu  = dAnbsduv[j][0];
        Real dAnbsdv  = dAnbsduv[j][1];
        Real dVnbsdu  = dVnbsduv[j][0];
        Real dVnbsdv  = dVnbsduv[j][1];

        // we can calculate the contact angle by finding the dgamma_gamma where dEdv is 0
        Real ca = 1.0 / dAnbsdv * (-dAdv + 2.0f * k0 * (dVdv + dVnbsdv));
        ca_list.push_back(ca);
    }
}

MeshTools::Real MeshTools::CalculateBoundaryAverageHeight(Mesh& m){
    std::vector<int> BoundaryIndices;
    MeshTools::CalculateBoundaryVerticesIndex(m, BoundaryIndices);

    Real z_height = 0.f;
    const auto& verts = m.getvertices();

    for (int i=0;i<BoundaryIndices.size();i++){
        int ind = BoundaryIndices[i];
        z_height += verts[ind].position_[2];
    }

    z_height = z_height / (Real)BoundaryIndices.size();

    return z_height;
}

void MeshTools::CalculateAV(Mesh& m, Real& A, Real& V){
    std::vector<Real> vecArea;
    std::vector<Real3> Normal;
    A =  CalculateArea(m, vecArea, Normal);
    V =  CalculateVolumeDivergenceTheorem(m, vecArea, Normal);
}

MeshTools::Real MeshTools::CalculateVnbsUnderneath(Mesh& m, AFP_shape* s, int projected_plane, int v_num, bool useNumerical)
{
    // for the shape, we calculate s --> obtain u and v 
    std::vector<int> BoundaryIndices;
    std::vector<Real> ulist, vlist;
    MeshTools::CalculateBoundaryVerticesIndex(m,  BoundaryIndices);

    // in this function, we also ordered the boundaryindices
    MeshTools::FindBoundaryUV(m, ulist, vlist, BoundaryIndices, s, true);

    return MeshTools::CalculateVnbsUnderneath(m, s, BoundaryIndices, ulist, vlist, projected_plane, v_num, useNumerical);
}



MeshTools::Real MeshTools::CalculateVnbsUnderneath(Mesh& m, AFP_shape* s, std::vector<int>& BoundaryIndices, std::vector<Real>& ulist, \
                                  std::vector<Real>& vlist, int projected_plane, int num_v, bool useNumerical)
{
    // We assume that the ulist is ordered
    Real V = 0.0;
    int numBoundary = ulist.size();
    const auto& verts = m.getvertices();
    for (int i=0;i<vlist.size();i++){
        Real vstep = (Constants::PI/2 - vlist[i]) / num_v;
        Real u = ulist[i];
        Real du = ulist[(i+1) % numBoundary] - ulist[i];
        if (du < 0){
            du += 2 * Constants::PI;
        }

        int ind = BoundaryIndices[i];

        #pragma omp parallel 
        {
            Real V_local= 0.0;
            #pragma omp for
            for (int j=0;j<num_v;j++){
                Real3 drdu, drdv, pos;
                Real v = vlist[i] + j * vstep;

                drdu = s->drdu(u,v, useNumerical);
                drdv = s->drdv(u,v, useNumerical);
                pos  = s->calculatePos(u,v);

                Real3 cp = LinAlg3x3::CrossProduct(drdu, drdv);
                Real A_projected = std::abs(cp[projected_plane]) * du * vstep;

                V_local += A_projected * pos[projected_plane];
            }

            #pragma omp critical
            {
                V += V_local;
            }
        }
    }

    return V;
}

MeshTools::Real MeshTools::CalculateMaxCurvature(Mesh& m, const std::vector<int>& BoundaryIndices){
    Real Pbar=0.0, Abar=0.0;
    const auto& vertices = m.getvertices();
    const auto& faces    = m.getFaces();

    for (int i=0;i<BoundaryIndices.size();i++){
        int this_index = BoundaryIndices[i];
        int next_index = BoundaryIndices[(i+1) % BoundaryIndices.size()];
        Real3 dist_vec;
        Real dist;

        m.getVertexDistance(vertices[this_index], vertices[next_index], dist_vec, dist);

        Pbar += std::sqrt(std::pow(dist_vec[0],2) + std::pow(dist_vec[1],2));
    }

    // calculate Abar
    std::vector<Real> vecArea;
    std::vector<Real3> Normal;
    MeshTools::CalculateArea(m, vecArea, Normal);

    #pragma omp parallel
    {
        Real A_local=0.0;
        #pragma omp for
        for (int i=0;i<vecArea.size();i++){
            A_local += vecArea[i] * Normal[i][2];
        }

        #pragma omp critical
        {
            Abar += A_local;
        }
    }
    

    return Pbar / (2 * Abar);
}

MeshTools::Real MeshTools::CalculateEta(Mesh& m, Real& Pbar, Real& Abar, Real3 Boundary_center){
    // calculate the vertex normal of the mesh
    m.CalcVertexNormals();

    // we calculate the arrange boundary indices
    std::vector<int> BoundaryIndices;
    MeshTools::CalculateBoundaryVerticesIndex(m,  BoundaryIndices, true, Boundary_center);

    // then we calculate Pbar 
    const auto& verts = m.getvertices();
    std::vector<Real3> tdr;
    Pbar = MeshTools::CalculatePbar(m, BoundaryIndices,tdr);
    Real eta=0.0;
    for (int i=0;i<BoundaryIndices.size();i++){
        int vert_ind = BoundaryIndices[i];
        Real3 mvec   = LinAlg3x3::CrossProduct({0,0,1}, tdr[i]);
        Real mdotN   = LinAlg3x3::DotProduct(-1.0 * mvec, verts[vert_ind].normals_);
        eta += mdotN;
    }

    eta = eta / Pbar;

    Abar = MeshTools::CalculateAbar(m);

    return eta;
}

MeshTools::Real MeshTools::CalculatePbar(Mesh& m, const std::vector<int>& BoundaryIndices, std::vector<Real3>& tdr){
    tdr.clear();

    const auto& vertices = m.getvertices();
    Real Pbar = 0.0;

    for (int i=0;i<BoundaryIndices.size();i++){
        int this_index = BoundaryIndices[i];
        int next_index = BoundaryIndices[(i+1) % BoundaryIndices.size()];
        Real3 dist_vec;
        Real dist;

        m.getVertexDistance(vertices[this_index], vertices[next_index], dist_vec, dist);
        tdr.push_back(dist_vec);

        Pbar += std::sqrt(std::pow(dist_vec[0],2) + std::pow(dist_vec[1],2));
    }

    return Pbar;
}

MeshTools::Real MeshTools::CalculateAbar(Mesh& m){
    std::vector<Real> A;
    std::vector<Real3> N;
    MeshTools::CalculateTriangleAreasAndFaceNormals(m, A, N);

    Real Abar=0.0;
    for (int i=0;i<A.size();i++){
        Abar += A[i] * N[i][2];
    }

    return Abar;
}
