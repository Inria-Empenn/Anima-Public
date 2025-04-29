#include "animaTCKReader.h"
#include <sstream>
#include <limits>
#include <algorithm> // Pour std::all_of

namespace anima
{
    // --- Constructeur ---

    TckReader::TckReader(const std::filesystem::path& filename) : m_filename(filename)
    {
        std::ifstream file(m_filename, std::ios::binary);
        if (!file.is_open()) {
            throw std::runtime_error("TckReader Error: Could not open file: " + m_filename.string());
        }

        // Configurer le flux pour lever des exceptions en cas d'erreur
        file.exceptions(std::ifstream::failbit | std::ifstream::badbit);

        try
        {
            parseHeader(file);
            readStreamlinesData(file);
        }
        catch (const std::ifstream::failure& e)
        {
            throw std::runtime_error("TckReader I/O Error: " + std::string(e.what()) + " in file " + m_filename.string());
        }
        catch (const std::runtime_error& e)
        {
            // Re-lancer les erreurs de parsing ou de formatage
            throw std::runtime_error("TckReader Error in file " + m_filename.string() + ": " + e.what());
        }
        // Le fichier est fermé automatiquement par le destructeur de std::ifstream (RAII)
    }

    // --- Accesseurs ---

    const std::map<std::string, std::string>& TckReader::getHeader() const noexcept
    {
        return m_header;
    }

    void TckReader::Update()
    {
        //TODO

    }

    /**
     * @brief Implémentation de la conversion en vtkPolyData.
     */
    vtkSmartPointer<vtkPolyData> TckReader::convertToVtkPolyData() const
    {
        // 1. Créer les objets VTK principaux
        vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
        vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();

        // --- Optionnel: Définir la précision des points VTK ---
        // Par défaut, vtkPoints utilise des floats (VTK_FLOAT).
        // Si les données TCK sont en double et que la précision est critique, décommentez :
        // if (m_dataType == SupportedDataType::FLOAT64_LE) {
        //     points->SetDataType(VTK_DOUBLE);
        //     std::cout << "Info: Setting VTK points data type to double." << std::endl;
        // } else {
        //     points->SetDataType(VTK_FLOAT); // Ou laisser le défaut
        // }
        // Remarque: Même si TCK est Float32LE, on stocke en `double` dans `m_streamlines`.
        // Laisser vtkPoints en float par défaut est souvent suffisant.
        // Si vous voulez *forcer* float même si m_streamlines est double:
        // points->SetDataType(VTK_FLOAT);


        // 2. Itérer sur chaque streamline chargée dans m_streamlines
        for (const Streamline& streamline : m_streamlines)
        {
            // Ignorer les streamlines vides (ne devrait pas arriver si la lecture est correcte)
            if (streamline.empty())
            {
                continue;
            }

            // 3. Créer une liste d'IDs VTK pour définir la connectivité de cette ligne (streamline)
            vtkSmartPointer<vtkIdList> pointIdList = vtkSmartPointer<vtkIdList>::New();
            pointIdList->SetNumberOfIds(streamline.size()); // Pré-allocation pour la performance

            // 4. Itérer sur chaque point de la streamline actuelle
            for (size_t i = 0; i < streamline.size(); ++i)
            {
                const Point& currentPoint = streamline[i]; // currentPoint est std::array<double, 3>

                // Insérer le point dans la structure vtkPoints.
                // InsertNextPoint retourne l'ID (index) du point nouvellement ajouté.
                // Note: vtkPoints->InsertNextPoint prend des doubles, la conversion
                // depuis notre std::array<double, 3> est directe. Si vtkPoints
                // est de type float, VTK gère la conversion double -> float.
                vtkIdType vtkPointId = points->InsertNextPoint(currentPoint[0], currentPoint[1], currentPoint[2]);

                // Ajouter l'ID du point VTK à la liste pour la polyline actuelle.
                pointIdList->SetId(i, vtkPointId);
            }

            // 5. Une fois tous les points d'une streamline traités,
            // insérer la liste d'IDs comme une nouvelle cellule (une ligne/polyline)
            // dans le vtkCellArray.
            lines->InsertNextCell(pointIdList);
        }

        // 6. Associer les points et les définitions de lignes (cellules) à l'objet vtkPolyData
        polyData->SetPoints(points);
        polyData->SetLines(lines);

        // L'appelant reçoit le vtkSmartPointer, qui gère la durée de vie de l'objet polyData.
        return polyData;
    }



    const std::vector<TckReader::Streamline>& TckReader::getStreamlines() const noexcept
    {
        return m_streamlines;
    }

    const std::filesystem::path& TckReader::getFilePath() const noexcept
    {
        return m_filename;
    }

    const std::string& TckReader::getHeaderValue(const std::string& key) const
    {
        try
        {
            return m_header.at(key); // at() lève std::out_of_range si la clé n'existe pas
        }
        catch (const std::out_of_range&)
        {
            throw std::out_of_range("TckReader Error: Header key '" + key + "' not found in file " + m_filename.string());
        }
    }


    // --- Méthodes Privées ---

    void TckReader::parseHeader(std::ifstream& file) {
        std::string line;

        // 1. Vérifier le magic string
        if (!std::getline(file, line) || line != "mrtrix tracks")
        {
            throw std::runtime_error("Invalid TCK magic string.");
        }

        // 2. Lire les paires clé-valeur
        bool foundDataType = false;
        while (std::getline(file, line) && line != "END")
        {
            std::size_t colon_pos = line.find(':');
            if (colon_pos == std::string::npos)
            {
                // Ignorer les lignes vides ou mal formées avant END ? Ou erreur ?
                // Pour être strict, considérons cela comme une erreur.
                if (!line.empty())
                {
                    throw std::runtime_error("Malformed header line: '" + line + "'");
                }
                continue; // Ignorer ligne vide
            }

            std::string key = line.substr(0, colon_pos);
            // Supprimer les espaces éventuels après la clé
            key.erase(key.find_last_not_of(" \t") + 1);

            std::string value = line.substr(colon_pos + 1);
            // Supprimer les espaces éventuels avant la valeur
            value.erase(0, value.find_first_not_of(" \t"));

            m_header[key] = value;

            if (key == "datatype")
            {
                foundDataType = true;
                // Normaliser la casse pour la comparaison
                std::transform(value.begin(), value.end(), value.begin(), ::tolower);
                if (value == "float32le")
                {
                    m_dataType = SupportedDataType::FLOAT32_LE;
                    m_coordSize = sizeof(float);
                }
                else if (value == "float64le")
                {
                    m_dataType = SupportedDataType::FLOAT64_LE;
                    m_coordSize = sizeof(double);
                }
                else
                {
                    m_dataType = SupportedDataType::UNSUPPORTED;
                    // Ne pas lancer d'erreur ici, mais vérifier après la boucle
                    // On pourrait supporter d'autres types (BE) ici avec byte-swapping
                }
            }
        }

        // 3. Vérifier si le marqueur END a été trouvé
        if (line != "END")
        {
            throw std::runtime_error("Missing 'END' marker in header.");
        }

        // 4. Vérifier si le datatype a été trouvé et est supporté
        if (!foundDataType)
        {
            throw std::runtime_error("Header missing mandatory 'datatype' field.");
        }
        if (m_dataType == SupportedDataType::UNSUPPORTED)
        {
            throw std::runtime_error("Unsupported datatype: '" + m_header["datatype"] + "'. Only Float32LE and Float64LE are currently supported.");
        }

        // 5. Stocker la position de fin de l'en-tête (début des données binaires)
        m_dataOffset = file.tellg();
        if (m_dataOffset < 0)
        {
            throw std::runtime_error("Could not determine data offset after header.");
        }

        // 6. Calculer la taille d'un point
        m_pointSize = 3 * m_coordSize;
    }


    // --- Lecture des données binaires ---

    void TckReader::readStreamlinesData(std::ifstream& file)
    {
        // Se positionner au début des données binaires
        file.seekg(m_dataOffset);
        if (!file)
        {
            throw std::runtime_error("Failed to seek to data offset.");
        }

        // Optimisation possible: Réserver de la mémoire si le champ 'count' est présent et fiable
        if (m_header.count("count"))
        {
            try
            {
                size_t estimated_count = std::stoull(m_header["count"]);
                // Réserver une capacité raisonnable, mais ne pas faire confiance aveuglément
                m_streamlines.reserve(std::min((size_t)1000000, estimated_count));
            }
            catch (...)
            {
                // Ignorer si la conversion échoue, on ne réserve juste pas.
            }
        }


        // Boucle principale de lecture des streamlines
        bool read_ok = true;
        while (file.peek() != EOF && read_ok)
        {
            Streamline currentStreamline;
            // Potentiellement réserver une taille moyenne pour une streamline
            // currentStreamline.reserve(100);

            if (m_dataType == SupportedDataType::FLOAT32_LE)
            {
                read_ok = readSingleStreamline<float>(file, currentStreamline);
            }
            else if (m_dataType == SupportedDataType::FLOAT64_LE)
            {
                read_ok = readSingleStreamline<double>(file, currentStreamline);
            }
            else
            {
                // Ne devrait jamais arriver ici grâce aux vérifications dans parseHeader
                throw std::logic_error("Internal error: trying to read unsupported data type.");
            }

            // Ajouter la streamline si elle n'est pas vide (ou si elle est valide)
            // Une streamline peut être vide si la fin est détectée immédiatement
            if (read_ok && !currentStreamline.empty())
            {
                try
                {
                    m_streamlines.push_back(std::move(currentStreamline));
                }
                catch (const std::bad_alloc& e)
                {
                    throw std::runtime_error("Memory allocation failed while storing streamlines. File might be too large.");
                }
            }
            else if (!read_ok && file.eof() && !currentStreamline.empty())
            {
                // Fichier terminé au milieu d'une streamline ? Erreur ou format spécial?
                // On pourrait choisir de l'ignorer ou de la garder. Pour l'instant, on ignore.
                std::cerr << "Warning: TCK file ended unexpectedly within a streamline in file " << m_filename.string() << std::endl;
            }
        }
        // On vide la mémoire non utilisée si on avait réservé trop
        m_streamlines.shrink_to_fit();
    }


    // --- Template pour lire une streamline ---

    template <typename T>
    bool TckReader::readSingleStreamline(std::ifstream& file, Streamline& currentStreamline)
    {
        std::array<T, 3> buffer; // Buffer pour lire un point (x, y, z)

        while (true)
        {
            // Lire les 3 coordonnées du point
            file.read(reinterpret_cast<char*>(buffer.data()), sizeof(buffer));

            // Vérifier si la lecture a réussi (et si on n'est pas à la fin du fichier)
            if (!file)
            {
                // Si on atteint EOF exactement après avoir lu un point complet, gcount() sera sizeof(buffer).
                // Si on atteint EOF *pendant* la lecture d'un point, gcount() sera < sizeof(buffer).
                // Si une autre erreur survient, failbit ou badbit sera levé (géré par l'exception globale).
                if (file.eof() && file.gcount() == 0)
                {
                    // Fin de fichier propre, atteinte avant de commencer un nouveau point.
                    return true; // Indique qu'on n'a pas pu lire *ce* point, fin normale.
                }
                else if (file.eof() && file.gcount() < sizeof(buffer))
                {
                    // Fin de fichier inattendue au milieu d'un point.
                    return false; // Indique une fin anormale.
                }
                else
                {
                    // Autre erreur de lecture (sera normalement attrapée par l'exception du flux)
                    return false; // Indique une erreur.
                }
            }


            // Conversion vers le type de stockage interne (double)
            // Pas strictement nécessaire si Point utilisait T, mais pour l'interface actuelle :
            Point currentPoint;
            currentPoint[0] = static_cast<double>(buffer[0]);
            currentPoint[1] = static_cast<double>(buffer[1]);
            currentPoint[2] = static_cast<double>(buffer[2]);


            // Vérifier si c'est la fin de la streamline (NaN ou Inf)
            if (isEndOfStreamline(buffer))
            {
                return true; // Fin normale de cette streamline
            }

            // Ajouter le point à la streamline actuelle
            try
            {
                currentStreamline.push_back(currentPoint);
            }
            catch (const std::bad_alloc& e)
            {
                throw std::runtime_error("Memory allocation failed while reading a streamline point.");
            }
        }
    }


    // --- Vérification de fin de streamline ---

    // Spécialisation pour float
    bool TckReader::isEndOfStreamline(const std::array<float, 3>& coords) noexcept
    {
        // Vérifie si toutes les coordonnées sont NaN ou toutes sont infinies
        return (std::isnan(coords[0]) && std::isnan(coords[1]) && std::isnan(coords[2])) ||
            (std::isinf(coords[0]) && std::isinf(coords[1]) && std::isinf(coords[2]));
        // Alternative C++11:
        // return std::all_of(coords.begin(), coords.end(), [](float v){ return std::isnan(v); }) ||
        //        std::all_of(coords.begin(), coords.end(), [](float v){ return std::isinf(v); });
    }

    // Spécialisation pour double
    bool TckReader::isEndOfStreamline(const std::array<double, 3>& coords) noexcept
    {
        // Vérifie si toutes les coordonnées sont NaN ou toutes sont infinies
        return (std::isnan(coords[0]) && std::isnan(coords[1]) && std::isnan(coords[2])) ||
            (std::isinf(coords[0]) && std::isinf(coords[1]) && std::isinf(coords[2]));
        // Alternative C++11:
        // return std::all_of(coords.begin(), coords.end(), [](double v){ return std::isnan(v); }) ||
        //        std::all_of(coords.begin(), coords.end(), [](double v){ return std::isinf(v); });
    }


}; // namespace